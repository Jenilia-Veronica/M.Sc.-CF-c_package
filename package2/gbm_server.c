/*
 * GBM MONTE CARLO OBSERVATORY - C engine + HTML dashboard
 * ---------------------------------------------------------
 * This C program does ALL the maths (GBM paths, percentiles, VaR, options,
 * convergence) and also acts as a tiny web server on http://localhost:8080.
 * The browser page (index.html) only draws what the C engine sends as JSON.
 *
 * Linux / macOS : gcc -O2 -o gbmserver gbm_server.c -lm && ./gbmserver
 * Windows (MinGW): gcc -O2 -o gbmserver.exe gbm_server.c -lm -lws2_32 && gbmserver.exe
 * Keep index.html in the same folder and run from that folder.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#ifdef _WIN32
#include <winsock2.h>
typedef SOCKET sock_t;
#define CLOSESOCK closesocket
#else
#include <signal.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
typedef int sock_t;
#define CLOSESOCK close
#endif

/* ---------- random numbers: xoshiro256** + Box-Muller ---------- */
static uint64_t st[4];
static int have_spare = 0;
static double spare;
static uint64_t splitmix(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static void rng_seed(uint64_t s) {
    for (int i = 0; i < 4; i++) st[i] = splitmix(&s);
    have_spare = 0;
}
static uint64_t rng_next(void) {
    uint64_t r = rotl(st[1] * 5, 7) * 9, t = st[1] << 17;
    st[2] ^= st[0]; st[3] ^= st[1]; st[1] ^= st[2]; st[0] ^= st[3];
    st[2] ^= t; st[3] = rotl(st[3], 45);
    return r;
}
static double uniform(void) { return ((rng_next() >> 11) + 1) * (1.0 / 9007199254740992.0); }
static double normal(void) {
    if (have_spare) { have_spare = 0; return spare; }
    double u1 = uniform(), u2 = uniform(), r = sqrt(-2.0 * log(u1));
    spare = r * sin(6.283185307179586 * u2);
    have_spare = 1;
    return r * cos(6.283185307179586 * u2);
}

/* ---------- helpers ---------- */
typedef struct { char *s; size_t n, cap; } Buf;
static void bp(Buf *b, const char *fmt, ...) {          /* printf into a growing buffer */
    char tmp[512];
    va_list ap;
    va_start(ap, fmt);
    int len = vsnprintf(tmp, sizeof tmp, fmt, ap);
    va_end(ap);
    if (b->n + len + 1 > b->cap) {
        b->cap = (b->n + len + 1) * 2;
        b->s = realloc(b->s, b->cap);
    }
    memcpy(b->s + b->n, tmp, len + 1);
    b->n += len;
}
static double ncdf(double x) { return 0.5 * erfc(-x / 1.4142135623730951); }
static int cmpd(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
static double qget(const char *q, const char *key, double def) {   /* read ?key=value */
    char pat[32];
    snprintf(pat, sizeof pat, "%s=", key);
    const char *p = strstr(q, pat);
    while (p && p != q && p[-1] != '?' && p[-1] != '&') p = strstr(p + 1, pat);
    return p ? atof(p + strlen(pat)) : def;
}
static double clampd(double v, double lo, double hi) { return v < lo ? lo : (v > hi ? hi : v); }

/* ---------- the simulation: everything the dashboard shows ---------- */
static void simulate(Buf *o, double mu, double sg, double T, int N, double K, uint64_t seed) {
    enum { n = 120, NPATH = 60, BINS = 44 };
    const double S0 = 100.0, R = 0.05;
    rng_seed(seed);
    double dt = T / n, drift = (mu - 0.5 * sg * sg) * dt, vol = sg * sqrt(dt);
    double *P = malloc(sizeof(double) * (size_t)N * (n + 1)), *c = malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {                              /* exact GBM paths */
        double s = S0;
        P[(size_t)i * (n + 1)] = s;
        for (int j = 1; j <= n; j++) { s *= exp(drift + vol * normal()); P[(size_t)i * (n + 1) + j] = s; }
    }
    bp(o, "{\"bands\":[");                                     /* percentile bands per time step */
    for (int j = 0; j <= n; j++) {
        for (int i = 0; i < N; i++) c[i] = P[(size_t)i * (n + 1) + j];
        qsort(c, N, sizeof(double), cmpd);
        bp(o, "%s[%.2f,%.2f,%.2f,%.2f,%.2f]", j ? "," : "", c[(int)(.05 * N)], c[(int)(.25 * N)],
           c[N / 2], c[(int)(.75 * N)], c[(int)(.95 * N)]);
    }
    bp(o, "],\"paths\":[");                                    /* a few full paths to draw */
    for (int k = 0; k < NPATH && k < N; k++) {
        bp(o, "%s[", k ? "," : "");
        for (int j = 0; j <= n; j++) bp(o, "%s%.2f", j ? "," : "", P[(size_t)k * (n + 1) + j]);
        bp(o, "]");
    }
    /* c is now the sorted terminal prices */
    double lo = c[(int)(.005 * N)], hi = c[(int)(.995 * N)], bw = (hi - lo) / BINS;
    long cnt[BINS] = {0}, up = 0;
    double sum = 0, sq = 0, tail = 0;
    int k5 = (int)(.05 * N);
    for (int i = 0; i < N; i++) {
        sum += c[i]; sq += c[i] * c[i];
        if (c[i] > S0) up++;
        if (i <= k5) tail += c[i];
        int b = (int)((c[i] - lo) / bw);
        if (c[i] >= lo && b >= 0 && b < BINS) cnt[b]++;
    }
    bp(o, "],\"hist\":{\"lo\":%.3f,\"hi\":%.3f,\"c\":[", lo, hi);
    for (int b = 0; b < BINS; b++) bp(o, "%s%ld", b ? "," : "", cnt[b]);
    double mean = sum / N, sd = sqrt(fmax(sq / N - mean * mean, 0));
    bp(o, "]},\"st\":{\"mean\":%.3f,\"sd\":%.3f,\"up\":%.4f,\"q05\":%.3f,\"es\":%.3f}",
       mean, sd, (double)up / N, c[k5], S0 - tail / (k5 + 1));

    // options: antithetic Monte Carlo vs Black-Scholes 
    double a = (R - 0.5 * sg * sg) * T, b = sg * sqrt(T), df = exp(-R * T), cs = 0, cs2 = 0, ps = 0;
    const int M = 30000;
    for (int i = 0; i < M; i++) {
        double z = normal(), x = S0 * exp(a + b * z), y = S0 * exp(a - b * z);
        double cc = df * (fmax(x - K, 0) + fmax(y - K, 0)) / 2;
        cs += cc; cs2 += cc * cc;
        ps += df * (fmax(K - x, 0) + fmax(K - y, 0)) / 2;
    }
    double cm = cs / M, se = sqrt(fmax(cs2 / M - cm * cm, 0) / M);
    double d1 = (log(S0 / K) + (R + 0.5 * sg * sg) * T) / b, d2 = d1 - b;
    double bc = S0 * ncdf(d1) - K * df * ncdf(d2), bpt = K * df * ncdf(-d2) - S0 * ncdf(-d1);
    bp(o, ",\"opt\":{\"c\":%.4f,\"se\":%.4f,\"p\":%.4f,\"bc\":%.4f,\"bp\":%.4f}", cm, se, ps / M, bc, bpt);

    /* convergence of the plain Monte Carlo call price */
    bp(o, ",\"conv\":[");
    double s1 = 0, s2 = 0, nx = 10;
    int first = 1;
    for (int i = 1; i <= 100000; i++) {
        double v = df * fmax(S0 * exp(a + b * normal()) - K, 0);
        s1 += v; s2 += v * v;
        if (i >= nx) {
            double m = s1 / i;
            bp(o, "%s[%d,%.4f,%.4f]", first ? "" : ",", i, m, sqrt(fmax(s2 / i - m * m, 0) / i));
            first = 0;
            nx = ceil(nx * 1.15);
        }
    }
    bp(o, "]");
    free(P); free(c);
}

/* ---------- tiny HTTP server ---------- */
static void send_all(sock_t c, const char *p, size_t len) {
    while (len > 0) {
        int r = send(c, p, (int)(len > 65536 ? 65536 : len), 0);
        if (r <= 0) return;
        p += r; len -= (size_t)r;
    }
}
static void respond(sock_t c, int code, const char *type, const char *body, size_t len) {
    char h[256];
    snprintf(h, sizeof h, "HTTP/1.1 %d OK\r\nContent-Type: %s\r\nContent-Length: %lu\r\n"
             "Connection: close\r\nCache-Control: no-store\r\n\r\n", code, type, (unsigned long)len);
    send_all(c, h, strlen(h));
    send_all(c, body, len);
}

int main(void) {
#ifdef _WIN32
    WSADATA w;
    WSAStartup(MAKEWORD(2, 2), &w);
#else
    signal(SIGPIPE, SIG_IGN);
#endif
    sock_t s = socket(AF_INET, SOCK_STREAM, 0);
    int one = 1, port = 8080;
    setsockopt(s, SOL_SOCKET, SO_REUSEADDR, (const char *)&one, sizeof one);
    struct sockaddr_in a;
    memset(&a, 0, sizeof a);
    a.sin_family = AF_INET;
    a.sin_addr.s_addr = htonl(INADDR_LOOPBACK);          /* only this computer can connect */
    for (; port < 8100; port++) {
        a.sin_port = htons((unsigned short)port);
        if (bind(s, (struct sockaddr *)&a, sizeof a) == 0) break;
    }
    if (port >= 8100 || listen(s, 8) != 0) { printf("Could not open a port.\n"); return 1; }
    printf("\n  GBM Monte Carlo Observatory - C engine running\n");
    printf("  Open  http://localhost:%d  in your browser (opening it for you...)\n", port);
    printf("  Press Ctrl+C to stop.\n\n");
    char cmd[128];
#ifdef _WIN32
    snprintf(cmd, sizeof cmd, "start http://localhost:%d", port);
#elif defined(__APPLE__)
    snprintf(cmd, sizeof cmd, "open http://localhost:%d", port);
#else
    snprintf(cmd, sizeof cmd, "xdg-open http://localhost:%d >/dev/null 2>&1", port);
#endif
    if (system(cmd) != 0) printf("  (could not auto-open: open the address above yourself)\n");

    for (;;) {
        sock_t c = accept(s, NULL, NULL);
        char req[4096], path[2048];
        int r = recv(c, req, sizeof req - 1, 0);
        if (r <= 0) { CLOSESOCK(c); continue; }
        req[r] = 0;
        if (sscanf(req, "GET %2047s", path) != 1) { CLOSESOCK(c); continue; }
        if (strncmp(path, "/api/sim", 8) == 0) {
            double mu = clampd(qget(path, "mu", 0.10), -0.5, 1.0), sg = clampd(qget(path, "sg", 0.20), 0.01, 2.0);
            double T = clampd(qget(path, "T", 1.0), 0.1, 10.0), K = clampd(qget(path, "K", 105), 1, 1000);
            int N = (int)clampd(qget(path, "N", 3000), 100, 20000);
            uint64_t seed = (uint64_t)qget(path, "seed", 1);
            clock_t t0 = clock();
            Buf o = {0};
            simulate(&o, mu, sg, T, N, K, seed);
            double ms = 1000.0 * (clock() - t0) / CLOCKS_PER_SEC;
            bp(&o, ",\"ms\":%.1f}", ms);
            respond(c, 200, "application/json", o.s, o.n);
            printf("  C engine: %d paths x 120 steps  mu=%.0f%% sigma=%.0f%%  ->  %.1f ms\n", N, mu * 100, sg * 100, ms);
            free(o.s);
        } else if (strcmp(path, "/") == 0 || strncmp(path, "/index", 6) == 0) {
            FILE *f = fopen("index.html", "rb");
            if (!f) {
                const char *m = "index.html not found. Run gbmserver from the folder that contains it.";
                respond(c, 404, "text/plain", m, strlen(m));
            } else {
                fseek(f, 0, SEEK_END);
                long len = ftell(f);
                fseek(f, 0, SEEK_SET);
                char *buf = malloc(len + 1);
                size_t got = fread(buf, 1, len, f);
                fclose(f);
                respond(c, 200, "text/html; charset=utf-8", buf, got);
                free(buf);
            }
        } else {
            respond(c, 404, "text/plain", "Not found", 9);
        }
        CLOSESOCK(c);
    }
}
