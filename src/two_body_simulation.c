#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double x;
    double y;
} Vector;

typedef enum {
    RUNGE_KUTTA,
    EULER,
} OdeMethod;

typedef struct {
    double mass;
    Vector position;
} Planetary;

typedef struct {
    Planetary a;
    Planetary b;
    Vector relative_velocity;
} TwoBodyModel;

typedef struct {
    double eccentricity;
    double mass_ratio;
    unsigned int period;
    double delta_time;
    OdeMethod method;
} TwoBodyController;

typedef void (*OdeFunction)(TwoBodyController *controller, TwoBodyModel *model);

double dot(Vector *lhs, Vector *rhs);
double *get_vector_component(Vector *v, size_t index);
Vector add(Vector *lhs, Vector *rhs);
Vector sub(Vector *lhs, Vector *rhs);
Vector mul_scalar(Vector *v, double s);
double length(Vector *v);
Vector unit(Vector *v);
Vector direction(Vector *lhs, Vector *rhs);

// x'' = -x ( 1 + q ) / r^3
// y'' = -y ( 1 + q ) / r^3
void runge_kutta(
    TwoBodyController *controller,
    TwoBodyModel *model
);
void euler(TwoBodyController *controller, TwoBodyModel *model);

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

    TwoBodyModel model;
    model.a.mass = 1.0;
    model.b.mass = controller.mass_ratio;

    model.a.position.x = model.b.mass / (model.a.mass + model.b.mass);
    model.a.position.y = 0;

    model.b.position.x = -model.a.mass / (model.a.mass + model.b.mass);
    model.b.position.y = 0;

    model.relative_velocity.x = 0;
    model.relative_velocity.y = sqrt((1 + controller.mass_ratio) * (1 + controller.eccentricity));

    OdeFunction ode = 0;
    while (!ode) {
        printf("0. Runge-Kutta\n1. Euler's Method\nEnter method: ");
        scanf("%u", &controller.method);
        if (controller.method == RUNGE_KUTTA) {
            ode = runge_kutta;
        } else if (controller.method == EULER) {
            ode = euler;
        } else {
            printf("Wrong input.\n");
        }
    }

    FILE *f = fopen("vertex.buf", "w");
    if (f) {
        for (double t = 0; t < controller.period; t += controller.delta_time) {
            fprintf(f,
                "(%.5lf, %.5lf)         ",
                model.a.position.x,
                model.a.position.y
            );
            fprintf(f,
                "(%.5lf, %.5lf)\n",
                model.b.position.x,
                model.b.position.y
            );

            ode(&controller, &model); 
        }
        fclose(f);
    }
    return 0;
}

// Passed the test
double dot(Vector *lhs, Vector *rhs) {
    double x = lhs->x * rhs->x;
    double y = lhs->y * rhs->y;

    return x + y;
}

double *get_vector_component(Vector *v, size_t index) {
    if (index == 0) {
        return &v->x;
    } else if (index == 1) {
        return &v->y;
    } else {
        exit(23);
    }
}

Vector add(Vector *lhs, Vector *rhs) {
    double x = lhs->x + rhs->x;
    double y = lhs->y + rhs->y;

    Vector result = { x, y };
    return result;
}

void add_assign(Vector *lhs, Vector *rhs) {
    lhs->x += rhs->x;
    lhs->y += rhs->y;
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

Vector mul_scalar(Vector *v, double s) {
    double x = v->x * s;
    double y = v->y * s;

    Vector result = { x, y };
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

void runge_kutta(
    TwoBodyController *controller,
    TwoBodyModel *model
) {
    // Slope of r vector
    Vector k_r[4];
    // Slope of v vector
    Vector k_v[4];

    // Position Vector: r_b - r_a
    Vector r;
    // Velocity Vector: v_b - v_a
    Vector v;

    for (size_t i = 0; i < 4; i++) {
        if (i == 0) {
            r = sub(&model->b.position, &model->a.position);
            v = model->relative_velocity;
        } else {
            double coef;
            if (i == 1 || i == 2) {
                coef = 0.5;
            } else {
                coef = 1.0;
            }

            Vector r0 = sub(&model->b.position, &model->a.position);

            v.x = k_v[i - 1].x * (controller->delta_time * coef) + model->relative_velocity.x;
            v.y = k_v[i - 1].y * (controller->delta_time * coef) + model->relative_velocity.y;

            r.x = v.x * (controller->delta_time * coef) + r0.x;
            r.y = v.y * (controller->delta_time * coef) + r0.y;
        }
        double r_magnitude = length(&r);

        for (size_t j = 0; j < 2; j++) {
            *get_vector_component(k_v + i, j) =
                -*get_vector_component(&r, j)
                * (1 + controller->mass_ratio)
                / pow(r_magnitude, 3.0);

            *get_vector_component(k_r + i, j) = *get_vector_component(&v, j);
        }
    }

    Vector r0 = sub(&model->b.position, &model->a.position);
    for (size_t i = 0; i < 2; i++) {
        *get_vector_component(&r0, i) +=
            (1.0 / 6.0 * controller->delta_time)
                * (*get_vector_component(k_r, i)
                + 2.0 * *get_vector_component(k_r + 1, i)
                + 2.0 * *get_vector_component(k_r + 2, i)
                + *get_vector_component(k_r + 3, i));

        *get_vector_component(&model->a.position, i) =
            -(model->b.mass / (model->a.mass + model->b.mass))
            * *get_vector_component(&r0, i);

        *get_vector_component(&model->b.position, i) =
            (model->a.mass / (model->a.mass + model->b.mass))
            * *get_vector_component(&r0, i);

        *get_vector_component(&model->relative_velocity, i) +=
            (1.0 / 6.0 * controller->delta_time)
                * (*get_vector_component(k_v, i)
                + 2.0 * *get_vector_component(k_v + 1, i)
                + 2.0 * *get_vector_component(k_v + 2, i)
                + *get_vector_component(k_v + 3, i));
    }
}

void euler(
    TwoBodyController *controller,
    TwoBodyModel *model
) {
    Vector r = sub(&model->b.position, &model->a.position);
    double r_magnitude = length(&r);

    Vector acceleration = mul_scalar(&r, -(1 + controller->mass_ratio) / pow(r_magnitude, 3.0));
    Vector speed = model->relative_velocity;

    Vector new_position = mul_scalar(&speed, controller->delta_time);
    // r += v * t
    add_assign(&r, &new_position); 

    Vector new_speed = mul_scalar(&acceleration, controller->delta_time);
    // v += a * t
    add_assign(&model->relative_velocity, &new_speed);

    model->a.position = mul_scalar(&r, -(model->b.mass / (model->a.mass + model->b.mass)));
    model->b.position = mul_scalar(&r, (model->a.mass / (model->a.mass + model->b.mass)));
}
