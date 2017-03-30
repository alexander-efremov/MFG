#ifndef MFG_UTILS_H
#define MFG_UTILS_H

#include <float.h>
#include <stdio.h>
#include <math.h>
#include "consts.h"

inline void print_matrix_to_file(int n, int m, double *data, const char *filename) {
    FILE *f = fopen(filename, "w");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int k = i * n + j;
            fprintf(f, "%20.14le ", data[k]);
        }
        fprintf(f, "\n ");
    }
    fclose(f);
}

inline void print_int(const char *str, int i) {
    printf("%s %d\n", str, i);
}

inline void print_int_double(const char *str, int i, double d) {
    printf("%s %d %f\n", str, i, d);
}

inline void print_double(const char *str, double d) {
    printf("%s %f\n", str, d);
}

inline void print_double(const char *str, double d1, double d2) {
    printf("%s %f %f\n", str, d1, d2);
}

inline void print_float(const char *str, float d) {
    printf("%s %f\n", str, d);
}

inline void print_double_exp(const char *str, double d) {
    printf("%s %e\n", str, d);
}

inline void print_float_exp(const char *str, float d) {
    printf("%s %e\n", str, d);
}

inline void print_matrix(double *a, int n, int m, int precision = 8) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int k = i * n + j;
            switch (precision) {
                case 1:
                    printf("%.1le ", a[k]);
                    break;
                case 2:
                    printf("%.2le ", a[k]);
                    break;
                case 3:
                    printf("%.3le ", a[k]);
                    break;
                case 4:
                    printf("%.4le ", a[k]);
                    break;
                case 5:
                    printf("%.5le ", a[k]);
                    break;
                case 6:
                    printf("%.6le ", a[k]);
                    break;
                case 7:
                    printf("%.7le ", a[k]);
                    break;
                case 8:
                    printf("%.8le ", a[k]);
                    break;
                default:
                    break;
            }
        }
        printf("\n");
    }
}

inline void print_matrix1(double *a, int n, int m, int precision = 8) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int k = i * n + j;
            switch (precision) {
                case 1:
                    printf("%.1f ", a[k]);
                    break;
                case 2:
                    printf("%.2f ", a[k]);
                    break;
                case 3:
                    printf("%.3f ", a[k]);
                    break;
                case 4:
                    printf("%.4f ", a[k]);
                    break;
                case 5:
                    printf("%.5f ", a[k]);
                    break;
                case 6:
                    printf("%.6f ", a[k]);
                    break;
                case 7:
                    printf("%.7f ", a[k]);
                    break;
                case 8:
                    printf("%.8f ", a[k]);
                    break;
                default:
                    break;
            }
        }
        printf("\n");
    }
}

inline void print_matrix(double *a, int n, int m, const char *text, int precision = 8) {
    printf("\n%s\n", text);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int k = i * n + j;
            switch (precision) {
                case 1:
                    printf("%.1le ", a[k]);
                    break;
                case 2:
                    printf("%.2le ", a[k]);
                    break;
                case 3:
                    printf("%.3le ", a[k]);
                    break;
                case 4:
                    printf("%.4le ", a[k]);
                    break;
                case 5:
                    printf("%.5le ", a[k]);
                    break;
                case 6:
                    printf("%.6le ", a[k]);
                    break;
                case 7:
                    printf("%.7le ", a[k]);
                    break;
                case 8:
                    printf("%.8le ", a[k]);
                    break;
                default:
                    break;
            }
        }
        printf("\n");
    }
}

inline void print_vector(double *a, int n, int precision = 8) {
    for (int k = 0; k < n; ++k) {
        switch (precision) {
            case 1:
                printf("%.1f ", a[k]);
                break;
            case 2:
                printf("%.2f ", a[k]);
                break;
            case 3:
                printf("%.3f ", a[k]);
                break;
            case 4:
                printf("%.4f ", a[k]);
                break;
            case 5:
                printf("%.5f ", a[k]);
                break;
            case 6:
                printf("%.6f ", a[k]);
                break;
            case 7:
                printf("%.7f ", a[k]);
                break;
            case 8:
                printf("%.8f ", a[k]);
                break;
            default:
                break;
        }
    }
}

inline void print_vector(int *a, int n) {
    for (int k = 0; k < n; ++k) {
        printf("%d ", a[k]);
    }
    printf("\n");
}

inline double get_l1_norm_vec(int x_len, int y_len, double *data) { // new
    double r = 0.;
    for (int i = 0; i < x_len; ++i)
        for (int j = 0; j < y_len; ++j)
            r += fabs(data[y_len * i + j]);
    return r / (x_len * y_len);
}

inline double get_l1_norm(double hx, double hy, int x_len, int y_len, double *data) { // old
    double r = 0.;
    for (int i = 0; i < x_len; ++i)
        for (int j = 0; j < y_len; ++j)
            r += fabs(data[y_len * i + j]);
    return r * hx * hy;
}

inline double get_l1_norm(double h, int x_len, double *data) {
    double r = 0.;
    r += data[0] * h / 2.;
    for (int i = 1; i < x_len - 1; ++i)
        r += fabs(data[i]) * h;
    r += data[x_len - 1] * h / 2.;
    return r;
}

inline double get_l1_norm_int_trapezoidal(double hx, double hy, int x_len, int y_len, double *data) {
    double r = 0.;
    for (int i = 0; i < x_len; ++i)
        for (int j = 0; j < y_len; ++j)
            r += 0.25 * (fabs(data[y_len * i + j]) +
                         fabs(data[y_len * (i + 1) + j]) +
                         fabs(data[y_len * i + j + 1]) +
                         fabs(data[y_len * (i + 1) + j + 1])) * hx * hy;
    return r;
}


inline double get_l_inf_norm(int x_len, int y_len, double *data) {
    double max = FLT_MIN;
    for (int i = 0; i < x_len; ++i)
        for (int j = 0; j < y_len; ++j)
            if (fabs(data[y_len * i + j]) > max)
                max = fabs(data[y_len * i + j]);
    return max;
}

inline bool is_empty_file(FILE *f) {
    long savedOffset = ftell(f);
    fseek(f, 0, SEEK_END);
    if (ftell(f) == 0) {
        return true;
    }
    fseek(f, savedOffset, SEEK_SET);
    return false;
}

inline void append_statistics(int ox_len, int oy_len, double tau, int iterCount, double err_l1_vec,
                              double err_l1_tr, double res_inf, double *extrem, double *extrem_err, int time_steps) {
    FILE *file;
    const char *filename = "/home/jane/ClionProjects/fem_circle/statistics.dat";
    file = fopen(filename, "a");
    if (file == NULL) {
        perror("Error opening file.");
        return;
    }
    if (is_empty_file(file)) {
        fprintf(file, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "OX", "OY", "TAU", "ITERCOUNT", "L1ERR-VEC",
                "L1ERR-TR", "MAXRESIDUAL", "MIN_RHO", "MAX_RHO", "MIN_ERR", "MAX_ERR", "TIMESTP");
    }

    fprintf(file, "%d\t%d\t%le\t%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%d\n", ox_len, oy_len, tau, iterCount,
            err_l1_vec, err_l1_tr, res_inf, extrem[0], extrem[1], extrem_err[0], extrem_err[1], time_steps);

    fclose(file);
}

inline void append_statistics_expl(int ox_len, int oy_len, double tau, double err_l1_vec,
                                   double err_l1_tr, double *extrem, double *extrem_err, int time_steps) {
    FILE *file;
    const char *filename = "/home/jane/ClionProjects/fem_circle/statistics.dat";
    file = fopen(filename, "a");
    if (file == NULL) {
        perror("Error opening file.");
        return;
    }
    if (is_empty_file(file)) {
        fprintf(file, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "OX", "OY", "TAU", "L1ERR-VEC",
                "L1ERR-TR", "MIN_RHO", "MAX_RHO", "MIN_ERR", "MAX_ERR", "TIMESTP");
    }

    fprintf(file, "%d\t%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%d\n", ox_len, oy_len, tau,
            err_l1_vec, err_l1_tr, extrem[0], extrem[1], extrem_err[0], extrem_err[1], time_steps);

    fclose(file);
}

inline double calc_array_sum(double *a, int ox_len, int oy_len, bool isAbs) {
    double res = 0;
    for (int i = 0; i < ox_len; i++) {
        for (int j = 0; j < oy_len; j++) {
            if (isAbs)
                res += fabs(a[i * oy_len + j]);
            else
                res += a[i * oy_len + j];
        }
    }
    return res;
}

inline double calc_array_sum(int *grid, double *a, int ox_len, int oy_len, bool isAbs) {
    double res = 0;
    for (int i = 0; i < ox_len; i++) {
        for (int j = 0; j < oy_len; j++) {
            int lev = grid[i * oy_len + j];

            if (lev >= 0) {
                if (isAbs)
                    res += fabs(a[i * oy_len + j]);
                else
                    res += a[i * oy_len + j];
            }
        }
    }
    return res;
}

inline double *calc_array_extrems(double *a, int ox_len, int oy_len) {
    double *res = new double[2];

    double maxRes = a[0];
    double minRes = a[0];
    for (int i = 0; i < ox_len; ++i) {
        for (int j = 0; j < oy_len; ++j) {
            double val = a[i * oy_len + j];
            if (val > maxRes) maxRes = val;
            if (val < minRes) minRes = val;
        }
    }

    res[0] = minRes;
    res[1] = maxRes;
    return res;
}

inline void print_surface(const char *filename, int ox_len, int oy_len,
                          double hx, double hy, int t, double a, double c, double x0, double y0,
                          double tau, double u, double v, double *data) {
    char name[650];
    sprintf(name, "%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f_a=%f_c=%f.dat",
            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v, a, c);
    FILE *file = fopen(name, "w");
    fprintf(file, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = 'x' 'y' %s\nZONE T='SubZone'",
            filename);
    fprintf(file, "\nI=%d J=%d K=%d ZONETYPE=Ordered", oy_len + 1, ox_len + 1, 1);
    fprintf(file, "\nDATAPACKING=POINT\nDT=(SINGLE SINGLE SINGLE)");
    for (int i = 0; i < ox_len + 1; i++)
        for (int j = 0; j < oy_len + 1; j++)
            fprintf(file, "\n%-30.20g  %-30.20g %-30.20g", i * hx, j * hy,
                    data[(oy_len + 1) * i + j]);

    fclose(file);
}

inline void print_line_along_x(const char *filename, int ox_len, int oy_len,
                               double hx, double hy, int t, double a, double c, double x0, double y0,
                               double tau, double u, double v, double *data, int fixed_y) {
    char name[650];
    sprintf(name, "line_by_x_%s_y_fix=%d_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f_a=%f_c=%f.dat",
            filename, fixed_y, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v, a, c);
    FILE *file = fopen(name, "w");
    fprintf(file, "TITLE = \"XY LINE\"\nVARIABLES = \"x\", \"rho\"\nZONE T=\"Only Zone\",");
    fprintf(file, " I=%d, F=POINT", ox_len + 1);
    for (int i = 0; i < ox_len + 1; i++)
        fprintf(file, "\n%-30.20g  %-30.20g", i * hx, data[(oy_len + 1) * i + fixed_y]);

    fclose(file);
}

inline void print_line_along_y(const char *filename, int ox_len, int oy_len,
                               double hx, double hy, int t, double a, double c, double x0, double y0,
                               double tau, double u, double v, double *data, int fixed_x) {
    char name[650];
    sprintf(name, "line_by_y_%s_x_fix=%d_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f_a=%f_c=%f.dat",
            filename, fixed_x, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v, a, c);
    FILE *file = fopen(name, "w");
    fprintf(file, "TITLE = \"XY LINE\"\nVARIABLES = \"y\", \"rho\"\nZONE T=\"Only Zone\",");
    fprintf(file, " I=%d, F=POINT", oy_len + 1);
    for (int j = 0; j < oy_len + 1; j++)
        fprintf(file, "\n%-30.20g  %-30.20g", j * hy, data[(oy_len + 1) * fixed_x + j]);

    fclose(file);
}

inline bool is_corner_node(int i, int j, int nx, int ny) {
    return (i == 0 && j == 0) || (i == 0 && j == ny - 1) || (i == nx - 1 && j == 0) || (i == nx - 1 && j == nx - 1);
}

inline bool is_border_node(int i, int j, int nx, int ny) {
    if (is_corner_node(i, j, nx, ny)) return false;
    return (i > 0 && j == 0) || (i == 0 && j > 0) || (i == nx - 1 && j > 0) || (i > 0 && j == nx - 1);
}

#endif //MFG_UTILS_H