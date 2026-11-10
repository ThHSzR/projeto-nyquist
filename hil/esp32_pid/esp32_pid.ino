#include "protocolo.hpp"

Pid pid;
char entrada[96];
std::size_t usados = 0;
bool excesso = false;

void setup() {
    Serial.begin(115200);
}

void loop() {
    while (Serial.available()) {
        const char c = static_cast<char>(Serial.read());
        if (c == '\n') {
            entrada[usados] = '\0';
            char resposta[96];
            if (excesso) {
                Serial.println("ERRO");
            } else {
                responder(pid, entrada, resposta, sizeof(resposta));
                Serial.println(resposta);
            }
            usados = 0;
            excesso = false;
        } else if (c != '\r') {
            if (usados < sizeof(entrada) - 1) {
                entrada[usados++] = c;
            } else {
                excesso = true;
            }
        }
    }
}
