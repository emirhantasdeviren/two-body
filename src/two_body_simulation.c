#include <stdio.h>
#include <math.h>

typedef struct {
    double x;
    double y;
} Vector;

typedef struct {
    double mass;
    Vector position;
    Vector speed;
} Planetary;

typedef struct {
    Planetary a;
    Planetary b;
} TwoBodyModel;

typedef struct {
    double eccentricity;
    double mass_ratio;
    unsigned int period;
    double delta_time;
    unsigned int method;
} TwoBodyController;

double dot(Vector *lhs, Vector *rhs);
double length(Vector *v);
Vector unit(Vector *v);
Vector sub(Vector *lhs, Vector *rhs);
Vector direction(Vector *lhs, Vector *rhs);

// x'' = -x ( 1 + q ) / r^3
// y'' = -y ( 1 + q ) / r^3
double runge_kutta(
    double (*fprime)(double, double),
    double x0,
    double y0,
    double dx,
    double x_end
);

int main(void) {
    TwoBodyController controller;

    printf("Enter period (T): ");
    scanf("%u", &controller.period);
    printf("Enter delta time (dt): ");
    scanf("%lf", &controller.delta_time);
    printf("Enter mass ratio: ");
    scanf("%lf", &controller.mass_ratio);
    printf("Enter eccentricity: ");
    scanf("%lf", &controller.eccentricity);
    printf("0. Runge-Kutta\n1. Euler's Method\nEnter method: ");
    scanf("%u", &controller.method);

    TwoBodyModel model;
    model.a.mass = 1.0;
    model.b.mass = controller.mass_ratio;

    model.a.position.x = model.a.mass / (model.a.mass + model.b.mass);
    model.a.position.y = 0;

    model.b.position.x = model.b.mass / (model.a.mass + model.b.mass);
    model.b.position.y = 0;

    model.a.speed.x = 0;
    model.a.speed.y = sqrt((1 + controller.mass_ratio) * (1 + controller.eccentricity));

    model.b.speed.x = 0;
    model.b.speed.y = 0;

    return 0;
}

// Passed the test
double dot(Vector *lhs, Vector *rhs) {
    double x = lhs->x * rhs->x;
    double y = lhs->y * rhs->y;

    return x + y;
}

// Passed the test
Vector sub(Vector *lhs, Vector *rhs) {
    double x = lhs->x - rhs->x;
    double y = lhs->y - rhs->y;

    Vector result = {
        x,
        y
    };

    return result;
}

// Passed the test
double length(Vector *v) {
    return sqrt(pow(v->x, 2) + pow(v->y, 2));
}

// Passed the test
Vector unit(Vector *v) {
    double len = length(v);

    double x = v->x / len;
    double y = v->y / len;

    Vector result = { x, y };

    return result;
}

Vector direction(Vector *lhs, Vector *rhs) {
    Vector a = sub(lhs, rhs);
    Vector result = unit(&a);

    return result;
}
