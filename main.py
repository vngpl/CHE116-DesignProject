import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# FILM HEAT TRANSF. COEFF.
h = 15 # W / (m^2 * C)

# ALUMINUM PROPERTIES
k = 204 # W / (m * C) @ 20 C
DENSITY_AL = 2707 # kg/m^3
AL_PRICE_PER_KG = 2.63 # $ / kg (based on MarketInsider)

# SURROUNDING AND SURFACE TEMPS
T_INFINITY = 30 # C
T_ZERO = 85 # C
THETA = (T_ZERO - T_INFINITY)

# PLATE DIMENSIONS
LENGTH = 0.100 # m
WIDTH = 0.050 # m
# ASSUMED THICKNESS
Z = 0.010 # m

AREA_RECT = LENGTH * WIDTH # m^2

# CREATE AN ARRAY OF POSSIBLE ARRANGEMENTS (from 1x1 to 40x40 square array)
boxArray = np.arange(1, 41, 1)
# SQRT OF SQUARE ARRAY GIVES NUM SQUARES (AKA BOXES)
numBoxes = np.square(boxArray)
# AREA OF EACH SQUARE GIVEN BY RECT AREA OVER NUM SQUARES
areaBoxes = AREA_RECT / numBoxes
# SQRT OF BOX AREA GIVES LENGTH OF ONE SIDE OF SQUARES
boxLength = np.sqrt(areaBoxes)

# SPACINGS FROM 75% TO 95% OF SQUARE SIDE LENGTH
percentS = np.arange(0.75, 1, 0.05)
# VARYING PIN LENGTHS FROM 1MM to 10MM
lengthArray = np.arange(0.001, 0.011, 0.001)

# COMBINE THE ARRAYS INTO A SINGLE DATA ARRAY
data = np.column_stack((boxArray, numBoxes, areaBoxes, boxLength))

# REPLICATE VALUES TO ACCOUNT FOR EACH SPACING CONFIGURATION
data = np.tile(data, (np.size(percentS), 1))
spacing_rep = np.repeat(percentS, np.size(boxArray))
data = np.column_stack((data, spacing_rep))

# REPLICATE VALUES TO ACCOUNT FOR EACH LENGTH CONFIGURATION
data = np.tile(data, (np.size(lengthArray), 1))
length_rep = np.repeat(lengthArray, np.size(spacing_rep))
data = np.column_stack((data, length_rep))

# CREATE A TABLE WITH NECESSARY PARAMETERS
df = pd.DataFrame(data=data, columns=["Square Array", "# Pins", "A (m^2)", "X (m)", "% Spacing", "L (m)"])

# CALCULATE FIN AND BARE PLATE DIMENSIONS USING SPECIFIED FORMULAS
df["D (m)"] = df["X (m)"] * df["% Spacing"]
df["P (m)"] = np.pi * df["D (m)"]
df["A_F (m^2)"] = np.pi * (df["D (m)"] / 2)**2
df["V (m^3)"] = df["A_F (m^2)"] * df["L (m)"]
df["S (m)"] = df["X (m)"] - df["D (m)"]
df["A_BARE (m^2)"] = AREA_RECT - df["A_F (m^2)"]
df["m"] = np.sqrt((h * df["P (m)"]) / (k * df["A_F (m^2)"]))
df["M"] = np.sqrt(h * df["P (m)"] * k * df["A_F (m^2)"]) * THETA
df["NUMER."] = np.sinh(df["m"] * df["L (m)"]) + (h / (df["m"] * k)) * np.cosh(df["m"] * df["L (m)"])
df["DENOM."] = np.cosh(df["m"] * df["L (m)"]) + (h / (df["m"] * k)) * np.sinh(df["m"] * df["L (m)"])

# CALCULATE FIN, BARE PLATE, AND TOTAL HEAT TRANSFER
df["Q_B"] = h * df["A_BARE (m^2)"] * THETA
df["Q_F"] = df["# Pins"] * df["M"] * (df["NUMER."] / df["DENOM."])
df["Q_T (W/m*C)"] = df["Q_B"] + df["Q_F"]

# CALCULATE FIN, BARE PLATE, AND TOTAL MASS
df["M_F"] = df["# Pins"] * (df["V (m^3)"] * DENSITY_AL)
M_B = (AREA_RECT * Z) * DENSITY_AL
df["M_T (kg)"] = M_B + df["M_F"]

# CALCULATE COST BASED ON TOTAL MASS OF CONFIGURATION
df["$"] = df["M_T (kg)"] * AL_PRICE_PER_KG

table = df[["Square Array", "# Pins", "% Spacing", "L (m)", "D (m)", "S (m)", "Q_B", "Q_F", "Q_T (W/m*C)", "M_T (kg)", "$"]]

# FIND OPTIMAL CONFIGURATION WHERE FIN HEAT TRANSFER >= 50
table_f = table.loc[(table["Q_F"] >= 50) & (table["Q_T (W/m*C)"] <= 62.5)]
table_f = table_f.sort_values(by=["$", "Q_F", "Q_T (W/m*C)"], ascending=[True, True, False])
optimal_config = table_f.iloc[0]

# GET OPTIMAL L, D, AND Q
optimal_L = optimal_config[3]
optimal_D = optimal_config[4]
optimal_Q = optimal_config[7]

# LOOP THROUGH VARYING LENGTH CONFIGS TO PLOT D VS Q WITH VARYING L
for length, group in table.groupby("L (m)"):
    plt.scatter(group["D (m)"], group["Q_F"], label=f"Length = {length:.3f} m")
# PLOT OPTIMAL POINT ON D VS Q PLOT
plt.scatter(optimal_D, optimal_Q, label="Optimal Configuration", c='k')
plt.title("Fin Heat Transfer vs Diameter")
plt.xlabel("Diameter (m)")panel

plt.ylabel("Q_fin (W/(m^2*C))")
plt.xlim()
plt.legend()
plt.grid()

# PLOT DIAMETER VS LENGTH VS HEAT TRANSFER
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(table["D (m)"], table["L (m)"], table["Q_F"])
ax.scatter(optimal_D, optimal_L, optimal_Q, label="Optimal Configuration", c='r')
ax.set_xlabel('Diameter (m)')
ax.set_ylabel('Length (m)')
ax.set_zlabel('Q_F (W/(m^2*C))')
ax.set_box_aspect(aspect=None, zoom=0.9)
ax.legend()

plt.title('Fin Heat Transfer, Diameter, and Length')
plt.show()
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Extracting relevant columns for plotting
x = filtered_table['% Spacing']
y = filtered_table['# Pins']
z = filtered_table['L (m)']

# Creating the scatter plot
sc = ax.scatter(x, y, z, c=filtered_table['Q_F'], cmap='viridis')

# Adding color bar
cbar = plt.colorbar(sc, ax=ax, pad=0.1)
cbar.set_label('Heat Transfer Q_F (W)')

highlight_point = [0.75, 1600, 0.009]
ax.scatter(*highlight_point, color='red', s=30, label='Optimized Value')
# Labeling the axes
ax.set_xlabel('Spacing percentage')
ax.set_ylabel('Number of Fins')
ax.set_zlabel('Length (m)')

# Setting the title
plt.title('Heat Transfer (Q_F) for Different Spacing, Number of Fins, and Length')

ax.legend()

plt.show()
plt.figure(figsize=(10,6))
plt.scatter(df["# Pins"], df["S (m)"])
highlighted_point = (1600,0.0004419417)
plt.scatter(*highlighted_point, color='red', s=100, edgecolor='black', marker='X', label='Optimized Value')
plt.xlabel('# of pins')
plt.ylabel('S(m)')
plt.title("# pins by spacing")
plt.grid(True)
plt.legend()
plt.show()
plt.figure(figsize=(10, 6))
for percent in percentS:
subset = table[table["% Spacing"] == percent]
plt.scatter(subset["$"], subset["Q_T (W/m*C)"], label=f'{percent * 100:.0f}% Spacing')

# Adding the specific configuration point from the screenshot
specific_config = {
"Square Array": 40,
"# Pins": 1600,
"% Spacing": 0.75,
"L (m)": 0.009,
"D (m)": 0.0013258252,
"S (m)": 0.0004419417,
"Q_B": 4.1238610193,
"Q_F": 50.9771769438,
"Q_T (W/m*C)": 55.1010379631,
"M_T (kg)": 0.1891662185,
"$": 0.4975071546
}
plt.scatter(specific_config["$"], specific_config["Q_F"], color='k', marker = 'x', s = 40,label='optimized value')

plt.xlabel('Cost ($)')
plt.ylabel('Heat Transfer (Q_F) (W/m*C)')
plt.title('Cost vs heat Transfer for Different Spacing Percentages')
plt.legend()
plt.grid(True)
plt.show()
