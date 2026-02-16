from burgers_solution import BurgersSolution


if __name__ == "__main__":
    solution = BurgersSolution("sample_000000")
    print(solution)
    print(solution.get_u(10, 0.02 + (1e-5 / 2)))