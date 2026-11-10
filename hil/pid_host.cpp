#include "esp32_pid/protocolo.hpp"
#include <iostream>
#include <string>

int main() {
    Pid pid;
    std::string linha;
    while (std::getline(std::cin, linha)) {
        char resposta[96];
        responder(pid, linha.c_str(), resposta, sizeof(resposta));
        std::cout << resposta << std::endl;
    }
}
