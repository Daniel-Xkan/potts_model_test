import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from matplotlib.colors import to_rgb

def plot_epistatic_plane(xy_base, de1, de2, de12):
# Adjust the figure size to make the canvas larger for annotations
    # Define the vertices of the first triangle
    dde = de12 - (de1 + de2)
    plane_color = 'darkgreen'
    vertices = [[0, 0, 0], [0, xy_base, de1], [xy_base, 0, de2]]

    # Create a 3D plot
    fig = plt.figure()
    fig.set_size_inches(10, 10)
    ax = fig.add_subplot(111, projection='3d')

    # Add the first triangle to the plot
    triangle = Poly3DCollection([vertices], alpha=0.5, color=plane_color)
    ax.add_collection3d(triangle)

    # Set the limits and labels
    ax.set_xlim([0, xy_base])
    ax.set_ylim([xy_base, 0])
    ax.set_zlim([-xy_base, 0])
    ax.xaxis.set_ticks_position('upper')
    ax.yaxis.set_ticks_position('upper')
    ax.view_init(elev=5, azim=240)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('ΔE', rotation=90)

    # Define the vertices of the second triangle
    vertices2 = [[xy_base, 0, de2], [0, xy_base, de1], [xy_base, xy_base, de12]]

    # Add the second triangle to the plot
    color = 'red' if dde < 0 else 'blue'
    triangle2 = Poly3DCollection([vertices2], alpha=0.5, color=color)
    ax.add_collection3d(triangle2)

    # Plot straight lines from each vertex of the triangles to the xy surface
    for vertex in vertices + vertices2:
        x, y, z = vertex
        ax.plot([x, x], [y, y], [z, 0], alpha=0.5, color='gray', linestyle='--')

    # Define the vertices of the square
    square_vertices = [[0, 0, 0], [0, xy_base, 0], [xy_base, xy_base, 0], [xy_base, 0, 0]]

    # Add the square to the plot
    square = Poly3DCollection([square_vertices], alpha=0.5, color='gray', edgecolor='gray')
    ax.add_collection3d(square)

    # Add a triangle with vertices (0, xy_base, de1), (xy_base, 0, de2), (xy_base, xy_base, de1 + de2)
    vertices3 = [[0, xy_base, de1], [xy_base, 0, de2], [xy_base, xy_base, de1 + de2]]
    triangle3 = Poly3DCollection([vertices3], alpha=0.5, color=plane_color)
    ax.add_collection3d(triangle3)

    # add dde lable:
    # Draw a solid line between the x-coordinates of the second and third triangle
    line_color = 'blue' if dde > 0 else 'red'
    ax.plot([vertices2[2][0], vertices3[2][0]],  # x-coordinates
            [vertices2[2][1], vertices3[2][1]],  # y-coordinates
            [vertices2[2][2], vertices3[2][2]],  # z-coordinates
            color=line_color, linewidth=2)

    # Label the ΔΔE value
    # Label the ΔE values at the vertices of the second triangle
    for i, vertex in enumerate(vertices2):
        x, y, z = vertex
        label_color = 'blue' if z > 0 else 'red'
        ax.text(x, y, z - 1,  # Adjusted Z position to be slightly above the vertex
            f"ΔE: {z:.2f}", color=label_color, fontsize=10, ha='center')
        
    label_position = [(vertices2[2][0] + vertices3[2][0]) / 2,
                      (vertices2[2][1] + vertices3[2][1]) / 2,
                      min(vertices2[2][2], vertices3[2][2]) - 2]  # Adjusted Z position to be slightly below the bottom of the line
    ax.text(label_position[0], label_position[1], label_position[2],
            f"ΔΔE: {dde:.2f}", color=line_color, fontsize=10, ha='center')
    
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    plt.show()

def plot_epistatic_bar(title, de_list):
    # Define bar width and positions
    bar_width = 0.2
    indices = np.arange(len(de_list))

    # Initialize the plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot bars for de1, de2, de12, and additive expectation
    for i, des in enumerate(de_list):
        name_pairs, de1, de2, de12 = des
        additive_de = de1 + de2
        dde = de12 - additive_de

        # Single mutant 1 (de1)
        ax.bar(indices[i] - 1.5 * bar_width, de1, bar_width, color='#f08080', label='ΔE Single 1' if i == 0 else "")

        # Single mutant 2 (de2)
        ax.bar(indices[i] - 0.5 * bar_width, de2, bar_width, color='#cf6161', label='ΔE Single 2' if i == 0 else "")

        # Double mutant (de12)
        # ax.bar(indices[i] + 0.5 * bar_width, de12, bar_width, color='blue', label='Double' if i == 0 else "")
        # Double mutant (de12) with RGB color 0000cd
        ax.bar(indices[i] + 0.5 * bar_width, de12, bar_width, color='#0000CD', label='ΔE Double' if i == 0 else "")
        # Additive expectation (de1 + de2) on the same x-axis as de12
        ax.bar(indices[i] + 0.5 * bar_width, additive_de, bar_width, color='none', edgecolor='black', linestyle='--', label='Additive' if i == 0 else "")


        # # Add red curved arrow for epistatic deviation
        # arrow_color = 'red' if dde < 0 else 'blue'
        # ax.annotate('', xy=(indices[i] + 0.5 * bar_width, de12), xytext=(indices[i] + 1.5 * bar_width, additive_de),
        #     arrowprops=dict(arrowstyle='->', color=arrow_color, lw=1.5))

    # Add horizontal line at y=0
    ax.axhline(0, color='black', linewidth=0.8)

    # Customize the plot
    ax.set_xticks(indices)
    ax.set_xticklabels([name_pairs for name_pairs, _, _, _ in de_list], rotation=45, ha='right')
    ax.set_ylabel('Fitness Change (ΔE)')
    ax.set_title(title)
    ax.legend()

    plt.tight_layout()
    plt.show()

def epistatic_pie (flip,non_flip, compensatory, resue, non_compensatory):
    # Prepare inner (3-way) and outer (2-way) data, normalize to ensure they sum to full circle
    inner_vals = np.array([compensatory, resue, non_compensatory], dtype=float)
    outer_vals = np.array([flip, non_flip], dtype=float)

    # Avoid division by zero by providing a minimal default if all values are zero
    if inner_vals.sum() == 0:
        inner_vals = np.array([1.0, 0.0, 0.0])
    if outer_vals.sum() == 0:
        outer_vals = np.array([1.0, 0.0])

    inner_fracs = inner_vals / inner_vals.sum()
    outer_fracs = outer_vals / outer_vals.sum()

    fig, ax = plt.subplots(figsize=(6, 6))

    # Color palettes for inner and outer rings
    inner_colors = ['#2e7d32', '#ffb300', '#c62828']   # compensatory, rescue, non_compensatory
    outer_colors = ['#1976d2', '#9e9e9e']              # flip, non_flip

    startangle = 90

    # Outer ring (flip / non_flip)
    # Outer ring (flip / non_flip) with explicit percentages (sum to 100%)
    outer_pct = outer_fracs * 100
    outer_labels = [
        f'Flip ({flip}) - {outer_pct[0]:.1f}%',
        f'Non-flip ({non_flip}) - {outer_pct[1]:.1f}%'
    ]
    ax.pie(
        outer_fracs,
        radius=1.3,
        startangle=startangle,
        colors=outer_colors,
        labels=outer_labels,
        labeldistance=1.05,
        wedgeprops=dict(width=0.3, edgecolor='w')
    )

    # Inner ring (compensatory / rescue / non_compensatory)
    ax.pie(
        inner_fracs,
        radius=1.0,
        startangle=startangle,
        colors=inner_colors,
        labels=[f'Compensatory ({compensatory})', f'Rescue ({resue})', f'Non-comp ({non_compensatory})'],
        labeldistance=0.75,
        wedgeprops=dict(width=0.4, edgecolor='w'),
        autopct=lambda pct: f'{pct:.1f}%'
    )

    ax.set(aspect='equal')
    # plt.title('Epistatic distribution (outer: flip/non-flip, inner: categories)')
    plt.show()

def epistatic_pie2(compensatory_flip, compensatory_nonflip,
                   noncompensatory_flip, noncompensatory_nonflip,
                   rescue_flip, rescue_nonflip,
                   rof_flip, rof_nonflip):
    # Assemble values and labels in the desired order (flip then non-flip for each category)
    values = np.array([
        compensatory_flip, compensatory_nonflip,
        noncompensatory_flip, noncompensatory_nonflip,
        rescue_flip, rescue_nonflip,
        rof_flip, rof_nonflip
    ], dtype=float)

    labels = [
        f'Compensatory (flip) ({int(compensatory_flip)})',
        f'Compensatory (non-flip) ({int(compensatory_nonflip)})',
        f'Non-comp (flip) ({int(noncompensatory_flip)})',
        f'Non-comp (non-flip) ({int(noncompensatory_nonflip)})',
        f'Rescue (flip) ({int(rescue_flip)})',
        f'Rescue (non-flip) ({int(rescue_nonflip)})',
        f'ROF (flip) ({int(rof_flip)})',
        f'ROF (non-flip) ({int(rof_nonflip)})'
    ]

    # Avoid empty pie by providing a tiny nonzero default
    if values.sum() == 0:
        values = np.ones_like(values)

    # Base colors for each category (compensatory, non-comp, rescue, rof)
    base_colors = ['#2e7d32', '#c62828', '#ffb300', '#6a1b9a']

    # Helper to darken/lighten a color (flip darker, non-flip lighter)
    def shade(hexcol, factor):
        rgb = np.array(to_rgb(hexcol)) * factor
        rgb = np.clip(rgb, 0, 1)
        return tuple(rgb)

    # Build colors: for each base color append flip (darker) then non-flip (lighter)
    colors = []
    for base in base_colors:
        colors.append(shade(base, 0.7))   # flip (darker)
        colors.append(shade(base, 1.15))  # non-flip (lighter, clipped)

    fig, ax = plt.subplots(figsize=(7, 7))
    startangle = 90

    wedges, texts, autotexts = ax.pie(
        values,
        labels=labels,
        colors=colors,
        startangle=startangle,
        autopct=lambda pct: f'{pct:.1f}%',
        wedgeprops=dict(edgecolor='w')
    )

    # Improve text contrast for small wedges
    for txt in autotexts:
        txt.set_color('white')
        txt.set_fontsize(8)

    ax.set(aspect='equal')
    plt.tight_layout()
    plt.show()