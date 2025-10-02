# 1D Burgers' Equation Solver - Simple Explanation

### What is Burgers' equation used for?

- It is a basic math model that shows how things like speed or waves change over time and space.
- This equation helps us understand how that bump moves forward and how it smooths out or spreads over time.

### Explanation of variables:

- **dt**: How much time passes between each update
- **dx**: Distance between each point along the line
- **u**: Speed at each point on the line at the current time
- **nt**: Number of time steps (how many times we update speed)
- **nx**: Number of points along the line where we measure speed
- **un**: Copy of the previous speed values used for calculations
- **nu**: Controls how much the bump smooths out (like friction or resistance)

### Output:

- The program saves the final speeds at each point to a file.
- You can use this data to make a graph showing how the bump moved and changed.
