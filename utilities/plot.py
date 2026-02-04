import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch

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

def epistatic_pie2(name, compensatory_flip, compensatory_nonflip,
                   noncompensatory_flip, noncompensatory_nonflip,
                   rescue_flip, rescue_nonflip,
                   gof_flip, gof_nonflip):
    import matplotlib.pyplot as plt
    # New ordering: put Non-compensatory first, then Compensatory, then Rescue, then GOF
    outer_vals = np.array([
        noncompensatory_flip, noncompensatory_nonflip,
        compensatory_flip, compensatory_nonflip,
        rescue_flip, rescue_nonflip,
        gof_flip, gof_nonflip
    ], dtype=float)

    inner_vals = np.array([
        noncompensatory_flip + noncompensatory_nonflip,
        compensatory_flip + compensatory_nonflip,
        rescue_flip + rescue_nonflip,
        gof_flip + gof_nonflip
    ], dtype=float)

    # Avoid zero-sum
    if inner_vals.sum() == 0:
        inner_vals = np.ones_like(inner_vals)
    if outer_vals.sum() == 0:
        outer_vals = np.ones_like(outer_vals)

    # Percentages for legends
    inner_pct = inner_vals / inner_vals.sum() * 100
    outer_pct = outer_vals / outer_vals.sum() * 100

    # Base colors for categories (Non-comp, Compensatory, Rescue, GOF) - swapped first two
    #original
    # base_colors = ['#c62828', '#2e7d32', '#ffb300', '#6a1b9a']
    #IBM
    base_colors = ['#e45704','#dc257f','#785ef0','#658fff']

    # Outer (segment) colors: for each category produce [flip (darker), nonflip (lighter)]
    outer_colors = []
    for c in base_colors:
        outer_colors.append('mistyrose')
        outer_colors.append('lightgray')

    fig, ax = plt.subplots(figsize=(7, 7))
    startangle = 90
    # place the provided name in the center of the pie chart
    ax.text(0, 0, name, ha='center', va='center', fontsize=14, fontweight='bold')
    
    # Draw outer ring (category totals) in new order: Non-comp, Compensatory, Rescue, GOF
    ax.pie(
        inner_vals,
        radius=1.0,
        startangle=startangle,
        colors=base_colors,
        # labels=[
        #     f'Non-comp ({int(inner_vals[0])})',
        #     f'Compensatory ({int(inner_vals[1])})',
        #     f'Rescue ({int(inner_vals[2])})',
        #     f'GOF ({int(inner_vals[3])})'
        # ],
        
        labeldistance=1.02,
        wedgeprops=dict(width=0.4, edgecolor='w')
    )

    # Draw inner ring (flip / non-flip segments) matching the new outer_vals order
    wedges, texts, autotexts = ax.pie(
        outer_vals,
        radius=1.4,
        startangle=startangle,
        colors=outer_colors,
        # labels=[
        #     f'Non-comp flip ({int(noncompensatory_flip)})', f'Non-comp nonflip ({int(noncompensatory_nonflip)})',
        #     f'Comp flip ({int(compensatory_flip)})', f'Comp nonflip ({int(compensatory_nonflip)})',
        #     f'Rescue flip ({int(rescue_flip)})', f'Rescue nonflip ({int(rescue_nonflip)})',
        #     f'GOF flip ({int(gof_flip)})', f'GOF nonflip ({int(gof_nonflip)})'
        # ],
        labeldistance=0.7,
        # autopct=lambda pct: f'{pct:.1f}%',
        autopct=lambda pct: f'',
        wedgeprops=dict(width=0.4, edgecolor='w')
    )

    # Improve inner text contrast
    for t in autotexts:
        t.set_color('white')
        t.set_fontsize(8)

    ax.set(aspect='equal')
    plt.tight_layout()
    plt.show()

    # --- New separate legend/percentage figure ---
    # category legend (outer ring) - order matches base_colors / inner_vals
    category_handles = [Patch(facecolor=base_colors[i], edgecolor='w') for i in range(len(base_colors))]
    category_labels = [
        f'Non-comp: {int(inner_vals[0])} ({inner_pct[0]:.1f}%)',
        f'Compensatory: {int(inner_vals[1])} ({inner_pct[1]:.1f}%)',
        f'Rescue: {int(inner_vals[2])} ({inner_pct[2]:.1f}%)',
        f'GOF: {int(inner_vals[3])} ({inner_pct[3]:.1f}%)'
    ]

    # segment legend (inner ring flip/non-flip) - order matches outer_vals
    segment_handles = [Patch(facecolor=outer_colors[i], edgecolor='w') for i in range(len(outer_colors))]
    segment_labels = [
        f'Non-comp flip: {int(noncompensatory_flip)} ({outer_pct[0]:.1f}%)',
        f'Non-comp nonflip: {int(noncompensatory_nonflip)} ({outer_pct[1]:.1f}%)',
        f'Comp flip: {int(compensatory_flip)} ({outer_pct[2]:.1f}%)',
        f'Comp nonflip: {int(compensatory_nonflip)} ({outer_pct[3]:.1f}%)',
        f'Rescue flip: {int(rescue_flip)} ({outer_pct[4]:.1f}%)',
        f'Rescue nonflip: {int(rescue_nonflip)} ({outer_pct[5]:.1f}%)',
        f'GOF flip: {int(gof_flip)} ({outer_pct[6]:.1f}%)',
        f'GOF nonflip: {int(gof_nonflip)} ({outer_pct[7]:.1f}%)'
    ]

    fig_leg = plt.figure(figsize=(6, 4))
    ax_leg = fig_leg.add_subplot(111)
    ax_leg.axis('off')

    # place two legends on the same axes
    cat_legend = ax_leg.legend(category_handles, category_labels, loc='upper center',
                               bbox_to_anchor=(0.5, 0.95), ncol=1, frameon=False, title='Categories')
    seg_legend = ax_leg.legend(segment_handles, segment_labels, loc='center',
                               bbox_to_anchor=(0.5, 0.35), ncol=1, frameon=False, title='Flip / Non-flip')

    # keep the first legend
    ax_leg.add_artist(cat_legend)

    plt.tight_layout()
    plt.show()
