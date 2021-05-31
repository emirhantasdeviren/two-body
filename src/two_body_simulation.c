#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double x;
    double y;
} Vector;

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
    unsigned int method;
} TwoBodyController;

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

    model.a.position.x = -model.b.mass / (model.a.mass + model.b.mass);
    model.a.position.y = 0;

    model.b.position.x = model.a.mass / (model.a.mass + model.b.mass);
    model.b.position.y = 0;

    model.relative_velocity.x = 0;
    model.relative_velocity.y = sqrt((1 + controller.mass_ratio) * (1 + controller.eccentricity));

    for (double t = 0; t < controller.period; t += controller.delta_time) {
        printf("t = %lf\n", t);
        printf(
            "A:\tPosition: %lf, %lf\n",
            model.a.position.x,
            model.a.position.y
        );

        printf(
            "B:\tPosition: %lf, %lf\n\n",
            model.b.position.x,
            model.b.position.y
        );
        runge_kutta(&controller, &model); 
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
    Vector k[4];

    Vector temp_position;

    for (size_t i = 0; i < 4; i++) {
        if (i == 0) {
            temp_position = sub(&model->b.position, &model->a.position);
        } else {
            double coef;
            if (i == 1 || i == 2) {
                coef = 0.5;
            } else {
                coef = 1.0;
            }

            Vector relative_position = sub(&model->b.position, &model->a.position);

            Vector velocity = {
                k[i - 1].x * (controller->delta_time * coef) + model->relative_velocity.x,
                k[i - 1].y * (controller->delta_time * coef) + model->relative_velocity.y
            };

            temp_position.x = velocity.x * (controller->delta_time * coef) + relative_position.x;
            temp_position.y = velocity.y * (controller->delta_time * coef) + relative_position.y;
        }
        double r_magnitude = length(&temp_position);

        for (size_t j = 0; j < 2; j++) {
            *get_vector_component(k + i, j) =
                -*get_vector_component(&temp_position, j)
                * (1 + controller->mass_ratio)
                / pow(r_magnitude, 3.0);
        }
    }

    Vector relative_position = sub(&model->b.position, &model->a.position);
    for (size_t i = 0; i < 2; i++) {
        *get_vector_component(&model->relative_velocity, i) +=
            (1.0 / 6.0 * controller->delta_time)
                * (*get_vector_component(k, i)
                + 2.0 * *get_vector_component(k + 1, i)
                + 2.0 * *get_vector_component(k + 2, i)
                + *get_vector_component(k + 3, i));

        *get_vector_component(&relative_position, i) +=
            *get_vector_component(&model->relative_velocity, i)
            * controller->delta_time;

        *get_vector_component(&model->a.position, i) =
            -(model->b.mass / (model->a.mass + model->b.mass))
            * *get_vector_component(&relative_position, i);

        *get_vector_component(&model->b.position, i) =
            (model->a.mass / (model->a.mass + model->b.mass))
            * *get_vector_component(&relative_position, i);
    }
}
