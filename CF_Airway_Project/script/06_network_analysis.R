# -------------------------------
# Load Libraries
# -------------------------------
library(GENIE3)
library(igraph)
library(ggraph)
library(ggplot2)

# -------------------------------
# Load GENIE3 weight matrix
# -------------------------------
weight_matrix <- readRDS("results/tables/GRN_importance_scores/GRN_weights.rds")

# -------------------------------
# Convert to edge list and filter top edges
# -------------------------------
edges <- GENIE3::getLinkList(weight_matrix)

# Keep top 500 edges for visualization
edges_top <- edges[1:500, ]

# Create graph
g <- graph_from_data_frame(edges_top[, c("regulatoryGene", "targetGene")], directed = TRUE)

# -------------------------------
# Identify hub TFs (out-degree centrality)
# -------------------------------
deg <- degree(g, mode = "out")
hub_TFs <- sort(deg, decreasing = TRUE)[1:10]
cat("Top hub TFs:\n")
print(hub_TFs)

# -------------------------------
# Create folders if not exist
# -------------------------------
dir.create("results/figures/GRN_networks", recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# Visualize full network (top 500 edges)
# -------------------------------
plot_path <- "results/figures/GRN_networks/GRN_top500_edges.png"
png(plot_path, width = 1200, height = 1000)
plot(g, vertex.size = 5, vertex.label.cex = 0.7, main = "Top 500 GRN edges")
dev.off()

# -------------------------------
# Highlight top 10 hub TFs in network
# -------------------------------
V(g)$color <- ifelse(V(g)$name %in% names(hub_TFs), "red", "skyblue")
V(g)$size <- ifelse(V(g)$name %in% names(hub_TFs), 8, 4)

highlight_path <- "results/figures/GRN_networks/GRN_top500_edges_hubTFs.png"
png(highlight_path, width = 1200, height = 1000)
plot(g, vertex.label.cex = 0.8, main = "GRN Network - Top Hub TFs Highlighted")
dev.off()

# -------------------------------
# Optional: save interactive network using ggraph
# -------------------------------
# ggraph(g, layout = "fr") +
#   geom_edge_link(alpha = 0.5) +
#   geom_node_point(aes(color = name %in% names(hub_TFs), size = degree(g, mode="out"))) +
#   geom_node_text(aes(label = name), repel = TRUE) +
#   theme_void()
# ggsave("results/figures/GRN_networks/GRN_interactive_network.png", width = 12, height = 10)