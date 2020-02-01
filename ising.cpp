#include <cmath>
#include <cstdio>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <vector>
#include <termios.h>
#include <unistd.h>

termios set_attributes(){
    termios conf, conf_orig;

    tcgetattr(0, &conf);
    conf_orig = conf;
    conf.c_lflag    &= ~(ICANON | ECHO);
    conf.c_cc[VMIN]  = 1;
    conf.c_cc[VTIME] = 0;
    
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &conf);
    printf("\033[?25l");

    return conf_orig;
}

void restore_attributes(termios conf){
    printf("\033[?25h\n");
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &conf);
}

class ising2d_simulator{
    std::mutex mtx;
    std::mt19937 mt;
    size_t width, height, rows;
    double J, h;
    double T_max = 4.0 / log(1.0 + sqrt(2.0));
    size_t T_level;
    std::vector<int> spin;

    void clear(){
        printf("\033[2J");
    }
    void putstr(int r, int c, const char *str){
        printf("\033[%d;%dH%s", r, c, str);
    }

    void draw_frame(){
        std::lock_guard<std::mutex> lock(mtx);
        const char wall[] = "\u2592\u2592"; // "▒▒"

        clear();
        for(size_t x=3;x<=width+2;x+=2) putstr(1, x, wall);
        for(size_t y=1;y<=rows+2;++y){ putstr(y, 1, wall); putstr(y, width + 3, wall); }
        for(size_t x=3;x<=width+2;x+=2) putstr(rows + 2, x, wall);
        putstr(1, width + 7, "高温");
        putstr(rows + 2, width + 7, "低温");
        putstr(rows + 3, 1, "[↑][↓]: 温度変更    [q]: 終了");
    }
    void draw_temperature(){
        std::lock_guard<std::mutex> lock(mtx);
        const char block[] = "|\u2591\u2591|"; // "|░░|"

        putstr(rows + 1, width + 7, block);
        for(size_t t=0;t<rows-1;++t) putstr(rows - t, width + 7, t < T_level ? block : "|  |");
        putstr(rows + 3, 1, "");
        fflush(stdout);
    }
    void draw_spins(){
        std::lock_guard<std::mutex> lock(mtx);

        for(size_t y=0;y<rows;++y){
            putstr(y + 2, 3, "");
            for(size_t x=0;x<width;++x){
                int uh, lh;
                uh = spin[(2 * y + 0) * width + x];
                lh = (2 * y + 1 < height) ? spin[(2 * y + 1) * width + x] : -1;
                if(uh == 1 && lh == 1) printf("\u2588"); // "█"
                if(uh != 1 && lh == 1) printf("\u2584"); // "▄"
                if(uh == 1 && lh != 1) printf("\u2580"); // "▀"
                if(uh != 1 && lh != 1) printf(" ");      // " "
            }
        }
        putstr(rows + 3, 1, "");
        fflush(stdin);
    }
public:
    ising2d_simulator(size_t w, size_t h): width(w), height(h), rows((h + 1) / 2), T_level(rows - 1), spin(w * h){
        std::random_device dev;
        std::uniform_int_distribution<> dist(0, 1);
        mt.seed(dev());

        for(size_t i=0;i<w*h;++i) spin[i] = dist(mt) * 2 - 1;
        set_param();
    }

    void set_param(double J_=1.0, double h_=0.0){
        J = J_; h = h_;
    }
    void increment_temperature(){
        if(T_level + 1 < rows){ ++T_level; draw_temperature(); }
    }
    void decrement_temperature(){
        if(T_level > 0){ --T_level; draw_temperature(); }
    }

    void show(){
        setvbuf(stdout, nullptr, _IOFBF, BUFSIZ);
        draw_frame();
        draw_temperature();
        draw_spins();
    }
    void update(size_t n=0){
        std::uniform_int_distribution<int> dist_idx(0, width * height - 1);
        std::uniform_real_distribution<>   dist_prob(0.0, 1.0);
        auto idx = [&](int x, int y){
            x = (x + width ) % width;
            y = (y + height) % height;
            return y * width + x;
        };
        double T = T_max / (rows - 1) * T_level;

        if(n == 0) n = width * height;
        for(size_t q=0;q<n;++q){
            int pos = dist_idx(mt);
            auto [y, x] = div(pos, width);
            double deltaE = (2*J * (spin[idx(x, y-1)] + spin[idx(x-1, y)] + spin[idx(x+1, y)] + spin[idx(x, y+1)]) - 2 * h) * spin[pos];
            double p = std::min(1.0, exp(-deltaE / T));

            if(dist_prob(mt) < p) spin[pos] *= -1;
        }
        draw_spins();
    }
};

int main(int argc, char** argv){
    size_t L = 32;
    double J = 1.0, h = 0.0;
    if(argc == 2){ L = atoll(argv[1]); }
    else if(argc == 4){ L = atoll(argv[1]); J = atof(argv[2]); h = atof(argv[3]); }
    
    auto attr = set_attributes();

    ising2d_simulator s{L, L};
    s.set_param(J, h);
    s.show();

    bool flg = true;
    std::thread thread_output([&](){
        while(flg){
            s.update();
            std::this_thread::sleep_for(std::chrono::milliseconds(30));
        }
    });
    std::thread thread_input([&](){
        int key, rb;
        int stat = 0;
        while(true){
            rb = read(STDIN_FILENO, &key, 1);
            if(rb == -1) break;
            if(key == 'q') break;

            if(stat == 0 && key == '\033') ++stat;
            else if(stat == 1 && key == '[') ++stat;
            else if(stat == 2){
                if(key == 'A') s.increment_temperature();
                if(key == 'B') s.decrement_temperature();
                stat = 0;
            }
            else stat = 0;
        }
        flg = false;
    });

    thread_output.join();
    thread_input.join();

    restore_attributes(attr);
    return 0;
}