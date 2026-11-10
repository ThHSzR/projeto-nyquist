#pragma once
#include "pid.hpp"
#include <cstdio>
#include <cstring>

// O mesmo parser e o mesmo controlador são usados no ESP32 e no executável de bancada.
inline void responder(Pid& pid, const char* linha, char* resposta, std::size_t tamanho) {
    if (std::strcmp(linha, "RESET") == 0) {
        pid = Pid{};
        std::snprintf(resposta, tamanho, "OK");
        return;
    }
    long indice;
    double erro;
    char extra;
    if (std::sscanf(linha, "%ld %lf %c", &indice, &erro, &extra) != 2
        || indice < 0 || !std::isfinite(erro)) {
        std::snprintf(resposta, tamanho, "ERRO");
        return;
    }
    std::snprintf(resposta, tamanho, "%ld %.15g", indice, pid.step(erro));
}
