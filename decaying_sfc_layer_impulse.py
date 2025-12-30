import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0, y0

# Parameters
lambda_val = 0.1  # decay coefficient
z0 = 1.0          # boundary position
z_max = 20.0      # upper limit for plotting

# Create z array
z = np.linspace(z0, z_max, 1000)

# Compute Bessel function arguments
arg = 2 * np.sqrt(lambda_val * z)
arg0 = 2 * np.sqrt(lambda_val * z0)

# Compute eigenfunction: Z(z) = J_0(...) - (J_0(...)/Y_0(...)) * Y_0(...)
# With boundary condition Z(z_0) = 0
j0_z = j0(arg)
y0_z = y0(arg)
j0_z0 = j0(arg0)
y0_z0 = y0(arg0)

Z = j0_z - (j0_z0 / y0_z0) * y0_z

# Normalize for better visualization
Z_normalized = Z / np.max(np.abs(Z))

# Create figure with multiple subplots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Eigenfunction vs z
ax = axes[0, 0]
ax.plot(z, Z_normalized, 'b-', linewidth=2, label='Z(z)')
ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax.axvline(x=z0, color='r', linestyle='--', alpha=0.5, label=f'z₀ = {z0}')
ax.plot(z0, 0, 'ro', markersize=8, label='Boundary Z(z₀)=0')
ax.set_xlabel('z', fontsize=12)
ax.set_ylabel('Z(z) (normalized)', fontsize=12)
ax.set_title('Spatial Eigenfunction Structure', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

# Plot 2: Bessel functions comparison
ax = axes[0, 1]
ax.plot(z, j0_z, 'b-', linewidth=2, label=f'$J_0(2\\sqrt{{\\lambda z}})$')
ax.plot(z, y0_z, 'r-', linewidth=2, label=f'$Y_0(2\\sqrt{{\\lambda z}})$')
ax.axvline(x=z0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('z', fontsize=12)
ax.set_ylabel('Bessel function value', fontsize=12)
ax.set_title('Component Bessel Functions', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

# Plot 3: Argument of Bessel functions
ax = axes[1, 0]
arg_plot = 2 * np.sqrt(lambda_val * z)
ax.plot(z, arg_plot, 'g-', linewidth=2.5)
ax.axvline(x=z0, color='r', linestyle='--', alpha=0.5, label=f'z₀')
ax.set_xlabel('z', fontsize=12)
ax.set_ylabel('$2\\sqrt{\\lambda z}$', fontsize=12)
ax.set_title('Bessel Function Argument', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

# Plot 4: Transient solution with spatial structure
ax = axes[1, 1]
times = [0, 0.5, 1.0, 2.0, 5.0]
colors = plt.cm.viridis(np.linspace(0, 1, len(times)))

for t, color in zip(times, colors):
    q_transient = Z_normalized * np.exp(-lambda_val * 1 * t)  # k=1 for visualization
    ax.plot(z, q_transient, linewidth=2.5, label=f't = {t}', color=color)

ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax.axvline(x=z0, color='r', linestyle='--', alpha=0.5)
ax.set_xlabel('z', fontsize=12)
ax.set_ylabel('$w(z,t)$ (normalized)', fontsize=12)
ax.set_title('Transient Decay: $w(z,t) = Z(z) e^{-\\lambda kt}$', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig('eigenfunction_structure.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"Parameters:")
print(f"  λ = {lambda_val}")
print(f"  z₀ = {z0}")
print(f"  Bessel argument at z₀: 2√(λz₀) = {arg0:.4f}")
print(f"\nBoundary condition check:")
print(f"  Z(z₀) = {Z[0]:.2e} ≈ 0 ✓")

