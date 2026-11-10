#pragma once
#include <cmath>

class Pid {
public:
    static constexpr double dt = 0.02;

    double step(double error) {
        // C(s) = sqrt(17)*(s+1)^2 / [s*(s+4)].
        constexpr double scale = 4.123105625617661;
        constexpr double kp = scale * 7.0 / 16.0;
        constexpr double ki = scale / 4.0;
        constexpr double kd = scale * 9.0 / 64.0;
        integral += ki * dt * error;
        derivative = (derivative + 4.0 * kd * (error - previous)) / (1.0 + 4.0 * dt);
        previous = error;
        return kp * error + integral + derivative;
    }

private:
    double integral = 0.0;
    double derivative = 0.0;
    double previous = 0.0;
};
