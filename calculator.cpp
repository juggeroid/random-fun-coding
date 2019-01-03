#include <iostream>
#include <string>
#include <cmath>
#include <regex>
#include <iomanip>
#include <string_view>

class Calculator {
private:

  struct operator_t { char subtoken, priority, associativity; };

  long double result;
  std::string expression;
  std::string::const_iterator current = std::cbegin(expression);
  std::sregex_iterator match;
  std::stack<std::pair<operator_t, decltype (result)>> mid;
  std::stack<std::string::value_type> parentheses;
  const inline static std::regex pattern {"([0-9]*[.])?[0-9]+", std::regex_constants::optimize };

  operator_t get_operator() {
    switch (*current) {
      case '+': ++current; return {'+', 1, 'L'};
      case '-': ++current; return {'-', 1, 'L'};
      case '*': ++current; return {'*', 2, 'L'};
      case '/': ++current; return {'/', 2, 'L'};
      case '^': ++current; return {'^', 3, 'R'};
      case '%': ++current; return {'%', 2, 'L'};
      case ')': if (parentheses.empty() || parentheses.top() != '(') throw std::runtime_error("syntax error");
                else       return {'N', 0, 'L'};
      case '\0':           return {'N', 0, 'L'};
      default: throw std::runtime_error("syntax error");
   }
  }

  auto get_token() { std::smatch m = *match++; std::advance(current, m.length()); return std::stold(m.str()); }

  decltype (result) perform_operation(operator_t const& opr, decltype (result) lhs, decltype (result) rhs) const {
   switch (opr.subtoken) {
     case '+':      return lhs + rhs;
     case '-':      return lhs - rhs;
     case '/': if (!static_cast<bool>(rhs)) throw std::runtime_error("division by zero");
               else return lhs / rhs;
     case '*':      return lhs * rhs;
     case '%': if (!static_cast<bool>(rhs)) throw std::runtime_error("modulo by zero");
               else return std::fmod (lhs, rhs);
     case '^':      return std::pow  (lhs, rhs);
     default:       return 0;
   }
  }

  decltype (result) get_value() {
    decltype (result) value = 0;
    switch (*current) { // GNU case range extension is used, so it may not work under MSVS.
      case '0' ... '9':    value =  get_token();                                 break;
      case '(': ++current; parentheses.push('('); value = evaluate(); ++current; break;
      case '+': ++current; value =  get_value();                                 break;
      case '-': ++current; value = -get_value();                                 break;
      default: throw std::runtime_error("syntax error");
    }
    return value;
  }

  decltype (result) evaluate() {
    mid.push ({ {'N', 0, 'L'}, 0 });
    auto value = get_value();
    while (!mid.empty()) {
      auto opr { get_operator() };
      while (opr.priority < mid.top().first.priority || (opr.priority == mid.top().first.priority && opr.associativity == 'L')) {
        if (mid.top().first.subtoken == 'N') { mid.pop(); return value; }
        value = perform_operation(mid.top().first, mid.top().second, value);
        mid.pop();
      }
      mid.push({ opr, value }); value = get_value();
    }
    return 0;
  }

  Calculator(std::string_view str): expression(str), match {std::begin(expression), std::end(expression), pattern} {
    expression.erase(std::remove_if(std::begin(expression), std::end(expression), ::isspace), std::end(expression));
    result = evaluate();
  }

public:
  static auto calc(std::string_view expr) { Calculator c {expr}; return c.result; }
};

int main() {
  for (;;) {
    try {
      std::cout << "<calc> ";
      std::string input;
      std::getline(std::cin, input);
      auto a = Calculator::calc(input);
      std::cout << std::setprecision(18) << a << "\n";
    } catch (std::exception const& exception) {
        std::cerr << exception.what() << "\n";
    }
  }
}
