import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import argparse

VALID_TYPES = set(['RNA', 'WGS', 'ChIP'])
VALID_MODES = set(['bamcoverage', 'bamcompare', 'mbs', 'computematrix'])
TYPEDIC = {
    'human': '--',
    'triticum': ':'
}
CDIC = {
    'macos24.3:arm64': '#1f77b4',
    'rhel8.8:x86_64': '#ff7f0e',
    'sles15.6:x86_64': '#2ca02c'
}


parser = argparse.ArgumentParser(description="Plot benchmark results.")
parser.add_argument("-i", type=str, required=True, help="Path to the input CSV file.")
parser.add_argument("-o", type=str, required=True, help="Path to the output plot file.")
args = parser.parse_args()

def validate_df(input):
    a = pd.read_csv(input)
    a['dt3'] = pd.to_timedelta(a['dt3'])
    a['dt4'] = pd.to_timedelta(a['dt4'])
    a['speedup'] = a.dt3.dt.total_seconds() / a.dt4.dt.total_seconds()
    # Valid types.
    invalid_values = set(a['type']) - VALID_TYPES
    invalid_values = [i for i in invalid_values if 'bins_' not in i and 'gtf_' not in i]
    # Assert and print violations if any
    if invalid_values:
        raise ValueError(f"Column contains invalid values: {invalid_values}")

    # Valid modes.
    for i in a['mode'].unique():
        assert i in VALID_MODES, f'mode {i} is not allowed.'
        
    # Valid platforms.
    valid_platforms = ['macos24.3:arm64','rhel8.8:x86_64','sles15.6:x86_64']
    for i in a['platform'].unique():
        assert i in valid_platforms, f'platform {i} is not allowed.'
    return a

def plot_df(a, of):
    fig, ax = plt.subplots(nrows=2, ncols=4, figsize=(12,6), tight_layout=True)

    rix = 0
    cix = 0
    cix2 = 0
    for mode in ['bamcoverage', 'bamcompare', 'mbs', 'computematrix']:
        if mode == 'bamcoverage':
            for mtype in a[a['mode'] == mode]['type'].unique():
                for ix,r in a[(a['mode'] == mode) & (a['type'] == mtype)].iterrows():
                    ax[rix, cix].plot(
                        ['v3','v4'],
                        [r['dt3'].total_seconds()/60, r['dt4'].total_seconds()/60],
                        c=CDIC[r['platform']],
                        ls=TYPEDIC[r['organism']]
                    )
                    ax[rix, cix].set_title(f"{mode} - {mtype}")
                cix += 1
        else:
            rix = 1
            if mode == 'mbs':
                for ix,r in a[(a['mode'] == mode)].iterrows():
                    ax[rix, cix2].plot(
                        ['v3','v4'],
                        [r['dt3'].total_seconds()/60, r['dt4'].total_seconds()/60],
                        c=CDIC[r['platform']],
                        ls=TYPEDIC[r['organism']]
                    )
                ax[rix, cix2].set_title(f"{mode}")
                cix2 += 1
            else:
                for ix,r in a[(a['mode'] == mode)].iterrows():
                    ax[rix, cix2].plot(
                        ['v3', 'v4'],
                        [r['dt3'].total_seconds()/60, r['dt4'].total_seconds()/60],
                        c=CDIC[r['platform']],
                        ls=TYPEDIC[r['organism']]
                    )
                    ax[rix, cix2].set_title(f"{mode}")
                cix2 += 1
            
    sns.barplot(
        data=a,
        x='mode',
        y='speedup',
        hue='platform',
        ax=ax[rix, cix],
        palette=CDIC
    )
    ax[rix, cix].axhline(y=1, color='red', linestyle='dashed')
    ax[rix, cix].set_xlabel("")
    ticks = ax[rix, cix].get_xticks()
    ax[rix, cix].set_xticks(ticks)
    ax[rix, cix].set_xticklabels(ax[rix, cix].get_xticklabels(), rotation=90)
    ax[rix, cix].legend_.remove()

    ax[0, cix].set_xticks([])
    ax[0, cix].set_yticks([])
    ax[0, cix].set_frame_on(False)
    ax[0,0].set_ylabel("Walltime (minutes)")
    ax[1,0].set_ylabel("Walltime (minutes)")


    # Custom legend for organism
    legend_labels = list(CDIC.keys())
    legend_colors = list(CDIC.values())
    legend_patches1 = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10)
                    for color in legend_colors
    ]

    # Custom legend for platform
    legend_labels2 = list(TYPEDIC.keys())
    legend_styles = list(TYPEDIC.values())
    legend_patches2 = [
        plt.Line2D([0], [0], color='gray', linestyle=style, linewidth=2)
        for style in legend_styles
    ]

    # Combine and add to one legend
    ax[0, cix].legend(legend_patches1 + legend_patches2, legend_labels + legend_labels2, loc='center')

    fig.savefig(of, dpi=300)

if __name__ == '__main__':
    df = validate_df(args.i)
    plot_df(df, args.o)