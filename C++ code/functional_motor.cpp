#include <gpiod.h>
#include <iostream>
#include <thread>
#include <chrono>

// compile using: g++ -o functional_motor functional_motor.cpp -lgpiod -lpthread

using namespace std::chrono_literals;

struct Stepper {
    int step_pin;
    int dir_pin;
    int steps;
    const char* name;
};

void stepper_test(gpiod_chip* chip, const Stepper& s) {
    gpiod_line* step = gpiod_chip_get_line(chip, s.step_pin);
    gpiod_line* dir  = gpiod_chip_get_line(chip, s.dir_pin);

    if (!step || !dir) {
        std::cerr << "Failed to get GPIO line for " << s.name << "\n";
        return;
    }

    if (gpiod_line_request_output(step, "stepper_test", 0) < 0 ||
        gpiod_line_request_output(dir,  "stepper_test", 0) < 0) {
        std::cerr << "Failed to request lines for " << s.name << "\n";
        return;
    }

    std::cout << "Testing " << s.name << "...\n";

    // Forward direction
    gpiod_line_set_value(dir, 1);
    for (int i = 0; i < s.steps; i++) {
        gpiod_line_set_value(step, 1);
        std::this_thread::sleep_for(1ms);
        gpiod_line_set_value(step, 0);
        std::this_thread::sleep_for(1ms);
    }

    std::this_thread::sleep_for(200ms);

    // Reverse direction
    gpiod_line_set_value(dir, 0);
    for (int i = 0; i < s.steps; i++) {
        gpiod_line_set_value(step, 1);
        std::this_thread::sleep_for(1ms);
        gpiod_line_set_value(step, 0);
        std::this_thread::sleep_for(1ms);
    }

    gpiod_line_release(step);
    gpiod_line_release(dir);
}

int main() {
    // Change this if GPIO controller is different
    gpiod_chip* chip = gpiod_chip_open_by_name("gpiochip0");
    if (!chip) {
        std::cerr << "Failed to open gpiochip0\n";
        return 1;
    }

    Stepper steppers[] = {
        {9,  6, 10, "Port Brake"},
        {5,  8, 20, "Starboard Brake"},
        {11, 15, 30, "CW Airbrake"},
        {14, 7, 40, "CCW Airbrake"}
    };

    for (auto& s : steppers) {
        stepper_test(chip, s);
        std::this_thread::sleep_for(500ms);
    }

    gpiod_chip_close(chip);
    return 0;
}
