#include <iostream>
#include <fstream>
#include <iomanip>
#include <unistd.h>
#include <fcntl.h>
#include <linux/i2c-dev.h>
#include <linux/spi/spidev.h>
#include <termios.h>
#include <sys/ioctl.h>
#include <cstring>
#include <vector>

// compile using: g++ -o PCB_Functional PCB_Functional.cpp

// ------------------- Utility: I2C helpers -------------------
int i2c_read_register(int fd, uint8_t reg, uint8_t &val) {
    if (write(fd, &reg, 1) != 1) return -1;
    if (read(fd, &val, 1) != 1) return -1;
    return 0;
}

int i2c_write_register(int fd, uint8_t reg, uint8_t val) {
    uint8_t buf[2] = {reg, val};
    return (write(fd, buf, 2) == 2) ? 0 : -1;
}

// ------------------- SPI helpers -------------------
int spi_transfer(int fd, uint8_t *tx, uint8_t *rx, size_t len) {
    struct spi_ioc_transfer tr{};
    tr.tx_buf = (unsigned long)tx;
    tr.rx_buf = (unsigned long)rx;
    tr.len = len;
    tr.speed_hz = 500000;
    tr.bits_per_word = 8;
    return ioctl(fd, SPI_IOC_MESSAGE(1), &tr);
}

// ------------------- UART helpers -------------------
int uart_open(const char *dev, int baud=B9600) {
    int fd = open(dev, O_RDWR | O_NOCTTY | O_SYNC);
    if (fd < 0) return -1;
    struct termios tty{};
    tcgetattr(fd, &tty);
    cfsetospeed(&tty, baud);
    cfsetispeed(&tty, baud);
    tty.c_cflag = (tty.c_cflag & ~CSIZE) | CS8;
    tty.c_iflag = 0;
    tty.c_oflag = 0;
    tty.c_lflag = 0;
    tcsetattr(fd, TCSANOW, &tty);
    return fd;
}

// ------------------- Main -------------------
int main() {
    std::ofstream log("sensor_log.txt");
    if (!log.is_open()) {
        std::cerr << "Failed to open log file\n";
        return 1;
    }

    // ==== SPI: ADXL345 (example: read device ID reg 0x00) ====
    int spi_fd = open("/dev/spidev0.0", O_RDWR); // could also be /dev/spidev1.0
    if (spi_fd >= 0) {
        uint8_t tx[2] = {0x80 | 0x00, 0x00}; // read DEVID (0x00)
        uint8_t rx[2] = {};
        if (spi_transfer(spi_fd, tx, rx, 2) >= 0) {
            log << "ADXL345 DEVID: 0x" << std::hex << int(rx[1]) << "\n";
        }
        close(spi_fd);
    } else log << "ADXL345 not found\n";

    // ==== I2C bus 0 (CHECK USING ls /dev/i2c-*) ====
    int i2c_fd = open("/dev/i2c-0", O_RDWR);

    if (i2c_fd >= 0) {
        uint8_t val;

        // Gyro IAM-20380 (0x69) WHO_AM_I at 0x75
        ioctl(i2c_fd, I2C_SLAVE, 0x69);
        if (i2c_read_register(i2c_fd, 0x75, val) == 0)
            log << "Gyro WHO_AM_I: 0x" << std::hex << int(val) << "\n";

        // Magnetometer QMC5883L (0x0D) ID at 0x0D/0x0A..0x0C
        ioctl(i2c_fd, I2C_SLAVE, 0x0D);
        if (i2c_read_register(i2c_fd, 0x0A, val) == 0)
            log << "QMC5883L ID: 0x" << std::hex << int(val) << "\n";

        // Barometer DPS310 (0x77) PROD_ID at 0x0D
        ioctl(i2c_fd, I2C_SLAVE, 0x77);
        if (i2c_read_register(i2c_fd, 0x0D, val) == 0)
            log << "DPS310 PROD_ID: 0x" << std::hex << int(val) << "\n";

        // TMP1075 (0x48–0x4B), read temperature register (0x00)
        for (int addr = 0x48; addr <= 0x4B; addr++) {
            ioctl(i2c_fd, I2C_SLAVE, addr);
            uint8_t buf[2];
            uint8_t reg = 0x00;
            if (write(i2c_fd, &reg, 1) == 1 && read(i2c_fd, buf, 2) == 2) {
                int16_t raw = (buf[0] << 8 | buf[1]) >> 4;
                double tempC = raw * 0.0625;
                log << "TMP1075@" << std::hex << addr << " Temp: " 
                    << std::dec << tempC << " C\n";
            }
        }

        close(i2c_fd);
    } else log << "I2C bus open failed\n";

    // ==== UART5 GPS ====
    int uart_fd = uart_open("/dev/ttyS5", B9600); // check actual location using "dmesg | grep tty"
    if (uart_fd >= 0) {
        char buf[128];
        int n = read(uart_fd, buf, sizeof(buf)-1);
        if (n > 0) {
            buf[n] = 0;
            log << "GPS raw: " << buf << "\n";
        }
        close(uart_fd);
    } else log << "GPS UART open failed\n";

    log.close();
    std::cout << "Check complete, see sensor_log.txt\n";
    return 0;
}
