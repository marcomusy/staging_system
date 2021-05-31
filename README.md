# Embryonic Mouse Staging System

Given an experimental set of 2D lines + time, compute the average shapes for each of
the time points and extrapolate the shapes at intermediate time points.

---

![](https://user-images.githubusercontent.com/32848391/120231578-f8f6de80-c251-11eb-81e2-ffa894fd859e.gif)


## Workflow


- Align each group using rigid procrustes alignment (PA). BUT: PA allows for improper rotations (subsample mirroring) so we keep the Left/Right samples separate this would also serve as control

![](https://user-images.githubusercontent.com/32848391/120230949-a668f280-c250-11eb-9211-288ce8910bcc.png)

- Apply Iterative Closest Point, using PA output as baseline

- Extract mean outlines:

![](https://user-images.githubusercontent.com/32848391/120231748-5db23900-c252-11eb-88ae-af189dd5e406.png)

- Refit the whole experimental sample


## Scripts:

- `python stages_builder.py`
  creates the averages as described above

- plot resulting growth of the limb for the 3 datasets as a function of time:

![](https://user-images.githubusercontent.com/32848391/120232038-ff398a80-c252-11eb-8a10-1a3e89bb66c7.png)

- `python stages_fitter.py`
  creates the animation at the header of this page.






