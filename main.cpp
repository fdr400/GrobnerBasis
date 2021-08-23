#include "include/GrobnerBasis.h"
#include "include/modular.h"

#include <boost/algorithm/string.hpp>
#include <boost/rational.hpp>
#include <gflags/gflags.h>
#include <iostream>
#include <new>
#include <string>

DEFINE_string(vars_order, "", "порядок перменных");
DEFINE_string(order, "lex", "упорядочение мономов");
DEFINE_string(mode, "minimum_basis", "режим работы");
DEFINE_int64(base, 2, "модуль в поле вычетов");
DEFINE_string(optimization, "default", "оптимизация алгоритма Бухбергера");

using namespace Grobner;

#if FIELD == 1
using Field = double;
#elif FIELD == 2
using Field = boost::rational<int64_t>;
#elif FIELD == 3
using Field = Modular;
#endif

#define Monom Monom<Field>
#define Order Order<Field>
#define Polynomial Polynomial<Field>
#define PolynomialsSet PolynomialsSet<Field>

void handle_flags() {
    // parse vars_order flag.
    if (FLAGS_vars_order != "") {
        std::vector<std::string> result;
        boost::split(result, FLAGS_vars_order, boost::is_any_of(","));
        std::vector<Monom::VariableNumberType> set_order;
        for (const auto& str : result) {
            set_order.emplace_back(std::stoi(str));
        }
        Monom::InitOrder(set_order);
    }
    // set up order and optimization.
    Order::InitOrder(FLAGS_order);
    PolynomialsSet::SetOptimization(FLAGS_optimization);
    // set up modulo for Z_p if it was given.
#if FIELD == 3
    Modular::SetModulo(FLAGS_base);
#endif
};
std::string gen_rules();

void expand_modes();
void elementary_reduction();
void reduction();
void s_polynomial();
void buchberger();
void autoreduction();
void minimum_basis();
void expand_optimizations();

int main(int argc, char* argv[]) {
    // set up gflags.
    gflags::SetUsageMessage(gen_rules());
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    handle_flags();

    // run program in required mode.
    if (FLAGS_mode == "expand_modes") {
        expand_modes();
    } else if (FLAGS_mode == "elementary_reduction") {
        elementary_reduction();
    } else if (FLAGS_mode == "reduction") {
        reduction();
    } else if (FLAGS_mode == "s_polynomial") {
        s_polynomial();
    } else if (FLAGS_mode == "buchberger") {
        buchberger();
    } else if (FLAGS_mode == "autoreduction") {
        autoreduction();
    } else if (FLAGS_mode == "minimum_basis") {
        minimum_basis();
    } else if (FLAGS_mode == "expand_optimizations") {
        expand_optimizations();
    } else {
        std::cout << "\tУказан несуществующий режим.\n";
        expand_modes();
    }

    gflags::ShutDownCommandLineFlags();
    return 0;
}

std::string gen_rules() {
    std::vector<std::string> rules = {
        "Многочлены вводятся в формате LATEX (без пробелов), переменные -- в "
        "формате x_{n}^{k} [сначала номер, потом степень] (например, "
        "x_1x_{17}+24x_{57}^3x_1^{23}+x_{71}^{23})",
        "Чтобы ввести систему многочленов, нужно сначала ввести количество "
        "многочленов в ней, потом ввести каждый многочлен через пробел",
        "Чтобы получить список команд, установите значение флага --mode равным "
        "expand_modes",
        "Доступные упорядочения мономов: lex, grlex, grevlex, invlex (за них "
        "отвечает флаг --order)",
        "Чтобы изменить стандартный порядок переменных, нужно присвоить флагу "
        "--vars_order список из номеров переменных через запятую (например, "
        "--vars_order=3,1,2)",
        "Чтобы выбрать поле, при компиляции присвойте ключу -DFIELD значение 1 "
        "(real), 2 (rational) или 3 (z_p, в этом случае ещё нужно "
        "инициализиовать флаг --base при запуске)",
        "Для выбора оптимизации алгоритма Бухбергера присвойте флагу "
        "--optimization одно из следующих значений: default, "
        "do_not_repeat_computation, skip_reciprocal, fastest_one",
        "Чтобы узнать, что означает каждая оптимизация, запустите программу в "
        "режиме --mode=expand_optimizations",
    };

    std::string res;
    const size_t max_len = 80;
    for (size_t i = 0; i < rules.size(); ++i) {
        std::string buff = "\n\t" + std::to_string(i + 1) + ") ";
        std::stringstream st(rules[i]);
        std::string curr;
        while (st >> curr) {
            if (buff.size() + curr.size() + 1 <= max_len) {
                buff += curr + " ";
            } else {
                res += buff + "\n\t";
                buff = curr + " ";
            }
        }
        res += buff;
        if (i + 1 == rules.size()) {
            res.back() = '.';
        } else {
            res.back() = ';';
        }
    }

    return res;
}

void expand_modes() {
    std::cout << "\tДоступные режимы:\n";
    std::cout << "\t [elementary_reduction]: производит редукцию одного "
                 "многочлена относительно другого.\n";
    std::cout << "\t [reduction]: производит редукцию многочлена, "
                 "относительно системы многочленов.\n";
    std::cout << "\t [s_polynomial]: находит S-полином двух многочленов.\n";
    std::cout << "\t [buchberger]: применяет алгоритм Бухбергера к системе "
                 "многочленов.\n";
    std::cout << "\t [autoreduction]: производит авторедукцию над системой "
                 "многочленов.\n";
    std::cout << "\t [minimum_basis]: ищет минимальный базис Грёбнера "
                 "идеала, порождённого заданной системой многочленов.\n";
}
void elementary_reduction() {
    std::cout << "Введите многочлен, который хотите отредуцировать:\n";
    Polynomial p1;
    std::cin >> p1;
    std::cout << "\nВведите многочлен, относительно которого хотите "
                 "отредуцировать:\n";
    Polynomial p2;
    std::cin >> p2;
    std::cout << "\n";

    auto [done, result] = p1.MakeElementaryReduction(p2);
    if (!done) {
        std::cout << "Многочлен нередуцируем относительного данного "
                     "многочлена.\n";
    } else {
        std::cout << "Результат: ";
        std::cout << result << "\n";
    }
}
void reduction() {
    std::cout << "Введите многочлен, который хотите отредуцировать:\n";
    Polynomial p;
    std::cin >> p;
    std::cout << "\nВведите систему многочленов, относительно которой "
                 "хотите отредуцировать:\n";
    PolynomialsSet ps;
    std::cin >> ps;
    std::cout << "\n";

    auto result = ps.MakeReduction(p);
    std::cout << "Результат: ";
    std::cout << result << "\n";
}
void s_polynomial() {
    std::cout << "Введите первый многочлен:\n";
    Polynomial p1;
    std::cin >> p1;
    std::cout << "\nВведите второй многочлен:\n";
    Polynomial p2;
    std::cin >> p2;
    std::cout << "\n";

    auto result = GetSPolynomial(p1, p2);
    std::cout << "Результат: ";
    std::cout << result << "\n";
}
void buchberger() {
    std::cout << "Введите систему многочленов, базис Грёбнера которой "
                 "хотите найти:\n";
    PolynomialsSet ps;
    std::cin >> ps;
    std::cout << "\n";

    ps.Buchberger();
    std::cout << "Результат: ";
    std::cout << ps << "\n";
}
void autoreduction() {
    std::cout << "Введите систему многочленов, над которой хотите "
                 "произвести процесс авторедукции:\n";
    PolynomialsSet ps;
    std::cin >> ps;
    std::cout << "\n";

    ps.AutoReduction();
    std::cout << "Результат: ";
    std::cout << ps << "\n";
}
void minimum_basis() {
    std::cout << "Введите систему многочленов, минимальный базис Грёбнера "
                 "которой "
                 "хотите найти:\n";
    PolynomialsSet ps;
    std::cin >> ps;
    std::cout << "\n";

    ps.BuildMinimumBasis();
    std::cout << "Результат: ";
    std::cout << ps << "\n";
}
void expand_optimizations() {
    std::cout << "\tДоступные оптимизации:\n";
    std::cout << "\t [default]: реализуется стандартный алгоритм Бухбергера.\n";
    std::cout << "\t [do_not_repeat_computation]: для каждой пары "
                 "многочленов S-полином считается только один раз.\n";
    std::cout << "\t [skip_reciprocal]: предыдущая оптимизация + если "
                 "лидеры двух многочленов взаимнопросто, то редукция их "
                 "S-полинома не производится.\n";
    std::cout << "\t [fastest_one]: предыдущие два + при подсчёте "
                 "S-полиномов используется "
                 "критерий для пространства сизигий системы многочленов.\n";
}