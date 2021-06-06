# Two-Body Numerical Simulation

A simulation to predict the motion of two planetary objects which are viewed as points.

## C
To compile in C using `gcc`

```$ gcc -Wall -o two_body_simulation src/two_body_simulation.c -lm```

We compiled in C17 standard but it also compiles in C99.

To run executable

```$ ./two_body_simulation```

You can view the results by opening `out_c.txt` with your favorite text editor.

## Python

To run with Python

```python src/two_body_simulation.py```

We run with Python 3.9.5 we don't know the minimum Python version to run our program.

You can view the results by opening `out_py.txt` with your favorite text editor.

## Output

Output file format is following

```(x_a, y_a), (x_b, y_b)```

Where each line represents the position at time.

---

Emirhan Taşdeviren - 270206027

Arda Metin Akdar - 270206010

Emre Üstere - 260201078
