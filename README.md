# Stochastic Processes

This repository contains MATLAB, Python, and Jupyter notebook resources for exploring and simulating **stochastic processes**.  
It includes theoretical studies, numerical simulations, and comparisons between stochastic integration methods for dynamical systems.

---

## 📁 Repository Structure and Contents

| File | Description |
|:-----|:------------|
| `CubicgeodesicMPC.ipynb` | Exploration of cubic geodesic paths versus Laplace assumptions in Model Predictive Control (MPC). |
| `CubicvsLA.ipynb` | Comparison between cubic paths and Laplace approximations. |
| `KramersSS.m` | MATLAB script for steady-state analysis of Kramers equation. |
| `Kramers_euler_maruyama.m` | MATLAB simulation of Kramers dynamics with Euler-Maruyama method. |
| `Kramersequation.ipynb` | Jupyter notebook exploring stochastic Kramers equation. |
| `Morris_Lecar_SS.ipynb` | Steady-state and stochastic analysis of Morris-Lecar neuron model. |
| `OUEulerMaruyama.py` | Python simulation of an Ornstein-Uhlenbeck process using Euler-Maruyama method. |
| `OUProcessSS.m` | MATLAB steady-state analysis of an Ornstein-Uhlenbeck process. |
| `OUdeterministicsimu.ipynb` | Deterministic simulation of Ornstein-Uhlenbeck dynamics. |
| `OUprocess.py` | Python module for Ornstein-Uhlenbeck process simulations. |
| `Pendulumwithfriction.ipynb` | Stochastic simulation of a pendulum system with friction. |
| `cubiclaplaceassum.m` | MATLAB function related to cubic approximation or Laplace assumptions. |
| `lorentzstochastic.ipynb` | Stochastic simulation of the Lorenz system. |
| `ornstein_uhlenbeck_euler_maruyama.m` | MATLAB implementation of Euler-Maruyama simulation for OU process. |
| `ouMilstein_py.ipynb` | Python notebook applying Milstein's method for stochastic differential equations. |
| `stochasticsimusOU.ipynb` | Collection of Ornstein-Uhlenbeck stochastic simulations. |
| `citation.cff` | Citation file for properly referencing this work. |

---

## 📚 Topics Covered

- **Ornstein-Uhlenbeck processes** (stochastic and deterministic)
- **Kramers equation** and simulations
- **Euler-Maruyama** and **Milstein methods** for SDEs
- **Stochastic modeling** of mechanical systems (e.g., pendulum with friction, Lorenz attractor)
- **Neuron model** dynamics with stochasticity (Morris-Lecar model)
- **Model Predictive Control** concepts linked with stochastic approximations

---

## 🛠 Requirements

- **MATLAB** R2020b or newer
- **Python 3.8+**
- Key Python libraries:
  - `numpy`
  - `matplotlib`
  - `scipy`
  - `sympy` (for some symbolic computations)
  - `jupyter` for notebooks

Install Python libraries with:
```bash
pip install numpy matplotlib scipy sympy jupyter
```

---

## 🚀 How to Use

1. **Clone the repository**:
   ```bash
   git clone https://github.com/YOUR_USERNAME/StochasticProcesses.git
   ```
2. **Open Jupyter Notebooks** for Python-based simulations:
   ```bash
   jupyter notebook
   ```
3. **Run MATLAB scripts** directly for theoretical and numerical analysis.

---

## ✨ Highlights

- Cross-discipline study of stochastic processes (physics, biology, control).
- Hands-on simulation of stochastic differential equations.
- Blending deterministic and stochastic dynamics in modeling.

---

## 📜 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file if available.

---

## 👨‍💻 Author

Developed by Adrian Guel.
