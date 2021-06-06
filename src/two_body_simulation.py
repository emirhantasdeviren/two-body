from dataclasses import dataclass
import math

@dataclass(repr = False)
class Vector:
    x: float
    y: float

    def __repr__(self):
        return "({:.5f}, {:.5f})".format(self.x, self.y)


    def __add__(self, other):
        return Vector(self.x + other.x, self.y + other.y)


    def __sub__(self, other):
        return Vector(self.x - other.x, self.y - other.y)


    def __mul__(self, value: float):
        return Vector(self.x * value, self.y * value)


    def __rmul__(self, value: float):
        return Vector(self.x * value, self.y * value)


    def __mul__(self, value: int):
        return Vector(self.x * value, self.y * value)


    def __rmul__(self, value: int):
        return Vector(self.x * value, self.y * value)


    def __truediv__(self, value: float):
        return Vector(self.x / value, self.y / value)


    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        return self


    def __neg__(self):
        return Vector(-self.x, -self.y)


    def length(self):
        return math.sqrt(self.x * self.x + self.y * self.y)

class Planetary:
    def __init__(self, mass: float, position: Vector):
        self.mass = mass
        self.position = position


class TwoBodyModel:
    def __init__(self, a: Planetary, b: Planetary, relative_velocity: Vector):
        self.a = a
        self.b = b
        self.relative_velocity = relative_velocity


class TwoBodyController:
    def __init__(self, period: int, delta_time: float, eccentricity: float, mass_ratio: float, method: int):
        self.period = period
        self.delta_time = delta_time
        self.eccentricity = eccentricity
        self.mass_ratio = mass_ratio
        self.method = method

    def __init__(self):
        self.period = int(input("Enter period (T): "))
        self.delta_time = float(input("Enter delta time (dt): "))
        self.mass_ratio = float(input("Enter mass_ratio (q): "))
        self.eccentricity = float(input("Enter eccentricity: "))
        while self.eccentricity < 0 or self.eccentricity > 1:
            self.eccentricity = float(input("\nWrong input. Eccentricity must be between 0 and 1\nEnter eccentricity: "))
        self.method = int(input("\n0. Runge-Kutta\n1. Euler's Method\nEnter method: "))
        while self.method != 0 and self.method != 1:
            self.method = int(input("\nWrong input. Method must be either 0 or 1.\n0. Runge-Kutta\n1. Euler's Method\nEnter method: "))



    def init_model(self):
        m_a = 1
        m_b = self.mass_ratio

        a = Planetary(m_a, Vector(m_b / (m_a + m_b), 0))
        b = Planetary(m_b, Vector(-m_a / (m_a + m_b), 0))

        relative_velocity = Vector(0, math.sqrt((1 + self.mass_ratio) * (1 + self.eccentricity)))

        return TwoBodyModel(a, b, relative_velocity)

    
    def euler(self, model: TwoBodyModel):
        r = model.b.position - model.a.position

        a = -r * (1 + self.mass_ratio) / math.pow(r.length(), 3)
        v = model.relative_velocity

        r += v * self.delta_time
        v += a * self.delta_time

        model.a.position = -(model.b.mass / (model.a.mass + model.b.mass)) * r
        model.b.position = (model.a.mass / (model.a.mass + model.b.mass)) * r


    def runge_kutta(self, model: TwoBodyModel):
        k_r = [Vector(0, 0), Vector(0, 0), Vector(0, 0), Vector(0, 0)]
        k_v = [Vector(0, 0), Vector(0, 0), Vector(0, 0), Vector(0, 0)]

        r = Vector(0, 0)
        v = Vector(0, 0)

        for i in range(4):
            if i == 0:
                r = model.b.position - model.a.position
                v = model.relative_velocity
            else:
                coef = 0.0
                if i == 1 or i == 2:
                    coef = 0.5
                else:
                    coef = 1.0

                r0 = model.b.position - model.a.position

                v = k_v[i - 1] * (self.delta_time * coef) + model.relative_velocity
                r = v * (self.delta_time * coef) + r0


            k_v[i] = -r * (1 + self.mass_ratio) / math.pow(r.length(), 3)
            k_r[i] = v

        r0 = model.b.position - model.a.position
        r0 +=  (1 / 6 * self.delta_time) * (k_r[0] + 2 * k_r[1] + 2 * k_r[2] + k_r[3])

        model.a.position = r0 * -(model.b.mass / (model.a.mass + model.b.mass))
        model.b.position = r0 * (model.a.mass / (model.a.mass + model.b.mass))

        model.relative_velocity +=  (1 / 6 * self.delta_time) * (k_v[0] + 2 * k_v[1] + 2 * k_v[2] + k_v[3])


def main():
    controller = TwoBodyController()
    model = controller.init_model()

    count = int(controller.period / controller.delta_time) + 1
    with open('out_py.txt', 'w') as f: 
        for _ in range(count):
            f.write(str(model.a.position))
            f.write(str(model.b.position))
            f.write('\n\n')
            if controller.method == 0:
                controller.runge_kutta(model)
            else:
                controller.euler(model)


if __name__ == "__main__":
    main()
