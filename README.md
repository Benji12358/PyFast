# PyFast

PyFast is Python tool environment developed to post-process MULTIFAST simulations. It was first developed in order to study transitional flows and to compute statistics inside a turbulent spot. Now, this tool is capable of:

- computing and plotting various statistics for both the velocity and the temperature,
- plotting different 2D contours,
- plotting 3D spots,
- computing and saving in files Lambda2 (with ParaView compatibility)

## Installation

To use PyFast, you just need to clone this repository.

## Usage

To use PyFast, you can check examples provided. All utilities functions are defined in the folder Utils/, while the main and the settings file must be at its root.
The settings file (in which you specify what you want to post-process) has to be defined for your own usage. To do so, there are examples for both a transitional turbulent case (`postprocessSettings_turb_transition.py`) and a theoretical fully turbulent case (`postprocessSettings_theo.py`).
Then, you need to write your own main file and function. Here again, there are examples corresponding to both cases, respectively (`main_turb_transition.py`) and (`main_theo.py`).

## Contributing

```bash
to be continued ...
```

## License

```bash
to be continued ...
```
