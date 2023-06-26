import hdbscan
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def pathogen_cluster(args):

    embedding_df = pd.read_csv(args.embedding, index_col=0)

    clustering_parameters = {
        **({"min_cluster_size": args.min_size} if args.min_size is not None else {}),
        **({"min_samples": args.min_samples} if args.min_samples is not None else {}),
        **({"cluster_selection_epsilon": args.distance_threshold} if args.distance_threshold is not None else {})
    }

    clusterer = hdbscan.HDBSCAN(**clustering_parameters)
    clusterer_default = hdbscan.HDBSCAN()

    clusterer.fit(embedding_df)
    clusterer_default.fit(embedding_df)
    embedding_df["label"] = clusterer.labels_.astype(str)
    embedding_df["label_default"] = clusterer_default.labels_.astype(str)

    if args.output_figure is not None:
        
        plot_data = {
            "x": embedding_df.to_numpy()[:, 0],
            "y": embedding_df.to_numpy()[:, 1],
        }

        plot_data["cluster"] = clusterer.labels_.astype(str)

        plot_df = pd.DataFrame(plot_data)
        ax = sns.scatterplot(
            data=plot_df,
            x="x",
            y="y",
            hue="cluster",
            alpha=0.5,
        )
        plt.savefig(args.output_figure)
        plt.close()
    
    if args.output_dataframe is not None:
        embedding_df.to_csv(args.output_dataframe, index_label="strain")
        