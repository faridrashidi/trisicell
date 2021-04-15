human = [
    ("1", 249698942, 125000000),
    ("2", 242508799, 93300000),
    ("3", 198450956, 91000000),
    ("4", 190424264, 50400000),
    ("5", 181630948, 48400000),
    ("6", 170805979, 61000000),
    ("7", 159345973, 59900000),
    ("8", 145138636, 45600000),
    ("9", 138688728, 49000000),
    ("10", 133797422, 40200000),
    ("11", 135186938, 53700000),
    ("12", 133275309, 35800000),
    ("13", 114364328, 17900000),
    ("14", 108136338, 17600000),
    ("15", 102439437, 19000000),
    ("16", 92211104, 36600000),
    ("17", 83836422, 24000000),
    ("18", 80373285, 17200000),
    ("19", 58617616, 26500000),
    ("20", 64444167, 27500000),
    ("21", 46709983, 13200000),
    ("22", 51857516, 14700000),
    ("X", 156040895, 60600000),
    ("Y", 57264655, 12500000),
]
mouse = [
    ("1", 196283350, 125000000),
    ("2", 182113224, 93300000),
    ("3", 160039680, 91000000),
    ("4", 157219549, 50400000),
    ("5", 153573022, 48400000),
    ("6", 149736546, 61000000),
    ("7", 145617427, 59900000),
    ("8", 129401213, 45600000),
    ("9", 124595110, 49000000),
    ("10", 130694993, 40200000),
    ("11", 122082543, 53700000),
    ("12", 120129022, 35800000),
    ("13", 120421639, 17900000),
    ("14", 124902244, 17600000),
    ("15", 104043685, 19000000),
    ("16", 98207768, 36600000),
    ("17", 94987271, 24000000),
    ("18", 90702639, 17200000),
    ("19", 61431566, 26500000),
    ("X", 171368232, 60600000),
    ("Y", 92500857, 12500000),
]


def _get_tree(
    newick,
    line_size,
    label_color,
    tiplab_size,
    inner_node_type,
    inner_node_size,
    distance_labels_to_bottom,
):
    cmd = f"""
    p <- ggtree(read.tree(text='{newick}'),
                layout='dendrogram',
                size={line_size}) %<+% info1 %<+% info2 +
    # geom_tippoint(aes(color={label_color}), size={tiplab_size}) +
    geom_tiplab(aes(color={label_color}),
                align=TRUE,
                linesize={line_size/4},
                angle=-90,
                hjust=-0.1,
                size={tiplab_size},
                fontface='bold') +
    scale_color_identity(guide='none') +
    geom_label2(aes(x=branch, subset=!isTip,
                    label={inner_node_type.lower()}_label),
                fill='white',
                color='black',
                # parse=TRUE,
                size={inner_node_size},
                label.padding=unit(0.1,'lines')) +
    xlim_tree({distance_labels_to_bottom})
    """
    return cmd


def _add_barplot(column_name, column_color, height=0.4):
    cmd = f"""
    p2 <- ggplot(data=info1) +
    geom_bar(aes(x=id_index, y={column_name}, fill={column_color}),
                 width=0.85, stat='identity', show.legend=FALSE) +
    scale_fill_gradient2(low='#E41A1C', midpoint=0, mid='white', high='#377EB8') +
    # scale_fill_identity(guide='none') +
    # scale_color_identity(guide='none') +
    theme_cowplot() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.title=element_text(size=6),
          legend.text=element_text(size=5),
          axis.title.y=element_text(size=6, face='bold'),
          axis.text.y=element_text(size=5),
          axis.line=element_line(size=0.2),
          axis.ticks=element_line(size=0.1),
          axis.ticks.length=unit(0.02, "in"),
    ) +
    labs(y='{column_name}', x='')
    p <- p %>% insert_bottom(p2, height={height})
    """
    return cmd


def _add_chromplot_helper(name, length):
    cmd = f"""
    p2 <- ggplot(data=info1) +
    geom_rect(aes(xmin=id_index-0.3, xmax=id_index+0.3,
                  ymax={length/1000000}, ymin=0, fill='white', color='black'),
              show.legend=FALSE) +
    scale_fill_identity(guide='none') +
    scale_color_identity(guide='none') +
    theme_cowplot() +
    theme(axis.title.x=element_blank(),
          axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.title=element_text(size=6),
          legend.text=element_text(size=5),
          axis.title.y=element_text(size=6, face='bold', angle=0, vjust=0.6, hjust=2),
          axis.text.y=element_blank(),
          axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.ticks.length=unit(0.02, "in"),
    ) +
    labs(y='chr{name}', x='')
    p2 <- p2 + geom_segment(x=4-0.3, y=2, xend=4+0.3, yend=2,
                            size=0.1, color='blue')
    p <- p %>% insert_bottom(p2, height={length/500000000})
    """
    return cmd


def _add_chromplot():
    cmd = ""
    for m in mouse:
        cmd += _add_chromplot_helper(m[0], m[1])
    return cmd
