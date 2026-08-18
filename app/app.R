library(shiny)
library(ggplot2)

dat <- readRDS("data/hmp_v35_five_habitats_40_each.rds")
counts <- dat$counts
metadata <- dat$metadata
stopifnot(identical(colnames(counts), rownames(metadata)))

habitat_levels <- c("Gut", "Nasal", "Oral", "Skin", "Vaginal")
metadata$HABITAT <- factor(metadata$HABITAT, levels = habitat_levels)
has_sex <- "SEX" %in% names(metadata)
if (!has_sex) metadata$SEX <- "Not available"
metadata$SEX <- factor(metadata$SEX)
sex_choices <- levels(metadata$SEX)

habitat_colors <- c(Gut = "#E76F51", Nasal = "#A5A51B", Oral = "#2BBF88",
  Skin = "#2EA8DF", Vaginal = "#D86BE8")

rank_tables <- list(Genus = counts)
if (!is.null(dat$family_counts)) rank_tables$Family <- dat$family_counts
if (!is.null(dat$phylum_counts)) rank_tables$Phylum <- dat$phylum_counts
rank_choices <- intersect(c("Phylum", "Family", "Genus"), names(rank_tables))

relative_abundance <- function(x) sweep(x, 2, colSums(x), "/")

diversity_values <- function(x) {
  p <- relative_abundance(x)
  data.frame(
    Richness = colSums(x > 0),
    Shannon = apply(p, 2, function(z) -sum(z[z > 0] * log(z[z > 0]))),
    Simpson = apply(p, 2, function(z) 1 - sum(z^2))
  )
}

distance_matrix <- function(x, method) {
  samples <- t(relative_abundance(x))
  if (method == "Bray-Curtis")
    return(as.dist(as.matrix(dist(samples, method = "manhattan")) / 2))
  binary <- 1 * (samples > 0)
  m <- as.matrix(dist(binary, method = "manhattan"))
  richness <- rowSums(binary)
  denominator <- outer(richness, richness, "+") + m
  jaccard <- ifelse(denominator == 0, 0, 2 * m / denominator)
  diag(jaccard) <- 0
  as.dist(jaccard)
}

model_ss <- function(a, group) {
  groups <- split(seq_along(group), factor(group))
  sum(vapply(groups, function(i) sum(a[i, i, drop = FALSE]) / length(i), numeric(1)))
}

permanova_stats <- function(distance, group, permutations = 999, seed = 407) {
  d <- as.matrix(distance); n <- nrow(d)
  h <- diag(n) - matrix(1 / n, n, n)
  a <- -0.5 * h %*% (d^2) %*% h
  ss_total <- sum(diag(a)); group <- droplevels(factor(group))
  df_model <- nlevels(group) - 1; df_resid <- n - nlevels(group)
  observed_ss <- model_ss(a, group)
  observed_f <- (observed_ss / df_model) / ((ss_total - observed_ss) / df_resid)
  set.seed(seed)
  permuted_f <- replicate(permutations, {
    ss <- model_ss(a, sample(group))
    (ss / df_model) / ((ss_total - ss) / df_resid)
  })
  data.frame(Term = "Group", Df = df_model, R2 = observed_ss / ss_total, F = observed_f,
    `Permutation p` = (1 + sum(permuted_f >= observed_f)) / (permutations + 1), check.names = FALSE)
}

dispersion_stats <- function(distance, group, permutations = 999, seed = 407) {
  n <- attr(distance, "Size")
  fit <- cmdscale(distance, k = min(n - 1, 40), eig = TRUE, add = TRUE)
  coordinates <- fit$points; group <- droplevels(factor(group))
  distances_to_centroid <- function(g) {
    g <- factor(g); out <- numeric(nrow(coordinates))
    for (level in levels(g)) {
      i <- which(g == level); center <- colMeans(coordinates[i, , drop = FALSE])
      centered <- sweep(coordinates[i, , drop = FALSE], 2, center)
      out[i] <- sqrt(rowSums(centered^2))
    }
    out
  }
  f_stat <- function(g) {
    z <- distances_to_centroid(g)
    unname(summary(aov(z ~ factor(g)))[[1]][1, "F value"])
  }
  observed_f <- f_stat(group); set.seed(seed)
  permuted_f <- replicate(permutations, f_stat(sample(group)))
  data.frame(Test = "Homogeneity of dispersion", F = observed_f,
    `Permutation p` = (1 + sum(permuted_f >= observed_f)) / (permutations + 1), check.names = FALSE)
}

pairwise_matrix_table <- function(result) {
  m <- result$p.value
  if (is.null(m)) return(data.frame())
  rows <- which(!is.na(m), arr.ind = TRUE)
  if (!nrow(rows)) return(data.frame())
  data.frame(Group1 = colnames(m)[rows[, "col"]], Group2 = rownames(m)[rows[, "row"]],
    `Adjusted p` = m[rows], check.names = FALSE)
}

pairwise_beta_tests <- function(distance, group, permutations = 199) {
  group <- droplevels(factor(group)); pairs <- combn(levels(group), 2, simplify = FALSE)
  do.call(rbind, lapply(seq_along(pairs), function(k) {
    pair <- pairs[[k]]; keep <- group %in% pair
    sub_distance <- as.dist(as.matrix(distance)[keep, keep]); sub_group <- droplevels(group[keep])
    perma <- permanova_stats(sub_distance, sub_group, permutations, 407 + k)
    disp <- dispersion_stats(sub_distance, sub_group, permutations, 807 + k)
    data.frame(Comparison = paste(pair, collapse = " vs "), PERMANOVA_R2 = perma$R2,
      PERMANOVA_p = perma[["Permutation p"]], Dispersion_p = disp[["Permutation p"]])
  }))
}

composition_summary <- function(rank_counts, sample_metadata, top_n) {
  rel <- relative_abundance(rank_counts); overall <- rowMeans(rel)
  top_taxa <- names(sort(overall, decreasing = TRUE))[seq_len(min(top_n, length(overall)))]
  result <- do.call(rbind, lapply(levels(droplevels(sample_metadata$HABITAT)), function(h) {
    i <- sample_metadata$HABITAT == h; means <- rowMeans(rel[, i, drop = FALSE]); shown <- means[top_taxa]
    data.frame(HABITAT = h, TAXON = c(names(shown), "Other"),
      ABUNDANCE = c(shown, max(0, 1 - sum(shown))), stringsAsFactors = FALSE)
  }))
  result$HABITAT <- factor(result$HABITAT, levels = habitat_levels)
  result$TAXON <- factor(result$TAXON, levels = c(top_taxa, "Other")); result
}

sample_composition <- function(rank_counts, sample_metadata, top_n) {
  rel <- relative_abundance(rank_counts)
  top_taxa <- names(sort(rowMeans(rel), decreasing = TRUE))[seq_len(min(top_n, nrow(rel)))]
  out <- do.call(rbind, lapply(seq_len(ncol(rel)), function(j) {
    shown <- rel[top_taxa, j]
    data.frame(SAMPLE = colnames(rel)[j], HABITAT = sample_metadata$HABITAT[j],
      TAXON = c(top_taxa, "Other"), ABUNDANCE = c(shown, max(0, 1 - sum(shown))))
  }))
  out$HABITAT <- factor(out$HABITAT, levels = habitat_levels)
  out$TAXON <- factor(out$TAXON, levels = c(top_taxa, "Other")); out
}

analysis_index <- function(mode, habitats, sexes, single_habitat = NULL) {
  if (mode == "Habitat") metadata$HABITAT %in% habitats & metadata$SEX %in% sexes
  else metadata$HABITAT == single_habitat & metadata$SEX %in% sexes
}

metric_text <- list(
  Richness = tags$span(strong("Observed richness"), " counts detected genera. Every genus has equal weight, and sequencing depth can influence the result."),
  Shannon = tags$span(strong("Shannon diversity"), " combines richness and evenness. Rare and moderately abundant genera influence the value; larger values indicate more diversity."),
  Simpson = tags$span(strong("Simpson diversity (1 − D)"), " emphasizes common and dominant genera. Values nearer 1 indicate more diversity; dominance pushes the value toward 0.")
)

ui <- navbarPage(
  title = "HMP Microbial Ecology Lab", id = "main_nav",
  header = tags$head(tags$link(rel = "stylesheet", href = "styles.css")),

  tabPanel("Start", value = "Start", div(class = "page-wrap",
    div(class = "hero", tags$span(class = "eyebrow", "BIOL 407 • HUMAN MICROBIOME PROJECT"),
      h1("Ask a question. Make a choice. Defend a conclusion."),
      p("Build a short ecological analysis using 200 samples from five human microbial habitats.")),
    fluidRow(column(4, div(class = "stat-card", h3("200"), p("samples"))),
      column(4, div(class = "stat-card", h3("5"), p("habitats"))),
      column(4, div(class = "stat-card", h3(nrow(counts)), p("retained genera")))),
    div(class = "lesson-card", h2("Guided analysis builder"),
      selectInput("lab_question", "1. Choose one ecological question", choices = c(
        "Does alpha diversity differ among groups?" = "alpha",
        "Does community composition differ among groups?" = "beta",
        "Which taxa characterize different habitats?" = "composition",
        "Does a particular taxon differ among habitats?" = "taxon")),
      textAreaInput("hypothesis", "2. State your hypothesis before seeing results", rows = 3,
        placeholder = "Example: Oral communities will have greater Shannon diversity than skin communities because..."),
      selectInput("planned_analysis", "3. Choose an analysis", choices = NULL),
      div(class = "callout", uiOutput("analysis_feedback")),
      actionButton("begin_analysis", "Begin analysis", class = "btn-primary"))
  )),

  tabPanel("Data", value = "Data", div(class = "page-wrap",
    h2("Meet the dataset"),
    p(class = "lead", "First inspect sampling depth. Each point is one sample; the box summarizes the distribution within a habitat."),
    fluidRow(column(4, div(class = "control-card",
      checkboxGroupInput("data_habitats", "Habitats", habitat_levels, habitat_levels),
      checkboxGroupInput("data_sexes", "Participant sex", sex_choices, sex_choices))),
      column(8, div(class = "plot-card", plotOutput("depth_plot", height = 380),
        div(class = "figure-note", strong("How to read this graph: "),
          "the y-axis is original sequencing depth on a logarithmic scale. Higher points contain more reads. A wide box indicates greater variation among samples.")))),
    div(class = "table-card", h3("Dataset summary"), tableOutput("sample_summary"),
      tags$dl(class = "definitions",
        tags$dt("Samples"), tags$dd("Number of retained samples in that habitat."),
        tags$dt("Median original reads"), tags$dd("Middle sequencing depth before taxonomic filtering."),
        tags$dt("Median retained reads"), tags$dd("Middle number of reads represented in the teaching table."),
        tags$dt("Median genera"), tags$dd("Middle observed genus richness across samples.")))
  )),

  tabPanel("Alpha diversity", value = "Alpha diversity", div(class = "page-wrap",
    h2("Question 1: Does within-sample diversity differ?"), uiOutput("hypothesis_banner"),
    fluidRow(column(3, div(class = "control-card",
      radioButtons("alpha_group", "Compare", c("Habitats" = "Habitat", "Sex within one habitat" = "Sex")),
      conditionalPanel("input.alpha_group == 'Habitat'",
        checkboxGroupInput("alpha_habitats", "Habitats", habitat_levels, habitat_levels),
        checkboxGroupInput("alpha_sexes", "Include sex", sex_choices, sex_choices)),
      conditionalPanel("input.alpha_group == 'Sex'",
        selectInput("alpha_one_habitat", "Habitat", setdiff(habitat_levels, "Vaginal"))),
      selectInput("alpha_metric", "Diversity metric", c("Observed richness" = "Richness",
        "Shannon diversity" = "Shannon", "Simpson diversity" = "Simpson")),
      radioButtons("alpha_test_choice", "Global test", c("Kruskal–Wallis" = "KW", "One-way ANOVA" = "ANOVA")),
      checkboxInput("alpha_pairwise", "Include pairwise comparisons", TRUE),
      actionButton("run_alpha", "Run alpha analysis", class = "btn-primary"))),
      column(9, div(class = "explanation-card", uiOutput("alpha_metric_explanation")),
        conditionalPanel("input.run_alpha == 0", div(class = "waiting-card", "Choose groups, metric, and test, then run the analysis.")),
        uiOutput("alpha_results_ui"))),
    div(class = "interpret-card", h3("Connect the result to your hypothesis"),
      textAreaInput("alpha_interpretation", NULL, rows = 3,
        placeholder = "Was your hypothesis supported? Cite the plot, effect pattern, and test result."))
  )),

  tabPanel("Beta diversity", value = "Beta diversity", div(class = "page-wrap",
    h2("Question 2: Does group membership structure community composition?"), uiOutput("hypothesis_banner_beta"),
    fluidRow(column(3, div(class = "control-card",
      radioButtons("beta_group", "Compare", c("Habitats" = "Habitat", "Sex within one habitat" = "Sex")),
      conditionalPanel("input.beta_group == 'Habitat'",
        checkboxGroupInput("beta_habitats", "Habitats", habitat_levels, habitat_levels),
        checkboxGroupInput("beta_sexes", "Include sex", sex_choices, sex_choices)),
      conditionalPanel("input.beta_group == 'Sex'",
        selectInput("beta_one_habitat", "Habitat", setdiff(habitat_levels, "Vaginal"))),
      radioButtons("distance_choice", "Community distance", c("Bray–Curtis" = "Bray-Curtis", "Jaccard" = "Jaccard")),
      checkboxInput("beta_pairwise", "Include pairwise habitat tests", FALSE),
      actionButton("run_beta", "Run community analysis", class = "btn-primary"))),
      column(9, div(class = "explanation-card", uiOutput("distance_explanation")),
        conditionalPanel("input.run_beta == 0", div(class = "waiting-card", "Choose groups and a distance measure, then run the analysis.")),
        uiOutput("beta_results_ui"))),
    div(class = "interpret-card", h3("Connect the result to your hypothesis"),
      textAreaInput("beta_interpretation", NULL, rows = 3,
        placeholder = "Describe separation, PERMANOVA R² and p, dispersion, and whether the hypothesis was supported."))
  )),

  tabPanel("Composition", value = "Composition", div(class = "page-wrap",
    h2("Question 3: Which taxa characterize each habitat?"), uiOutput("hypothesis_banner_composition"),
    fluidRow(column(3, div(class = "control-card",
      radioButtons("tax_rank", "Taxonomic rank", choices = rank_choices),
      sliderInput("top_n", "Number of named taxa", 4, 15, 8),
      checkboxGroupInput("composition_habitats", "Habitats", habitat_levels, habitat_levels),
      checkboxGroupInput("composition_sexes", "Include sex", sex_choices, sex_choices),
      radioButtons("composition_view", "Display", c("Habitat means" = "Mean", "Individual samples" = "Samples")),
      actionButton("run_composition", "Build composition graph", class = "btn-primary"))),
      column(9, div(class = "explanation-card", strong("Relative abundance"),
        " is the proportion of assigned reads belonging to each taxon. Bars sum to 100%, so an increase in one group necessarily changes the proportions of others."),
        conditionalPanel("input.run_composition == 0", div(class = "waiting-card", "Choose rank, habitats, and display, then build the graph.")),
        uiOutput("composition_results_ui"))),
    div(class = "interpret-card", h3("Describe the ecological pattern"),
      textAreaInput("composition_interpretation", NULL, rows = 3,
        placeholder = "Which taxa dominate each habitat? What habitat conditions might explain the pattern?"))
  )),

  tabPanel("Compare a taxon", value = "Compare a taxon", div(class = "page-wrap",
    h2("Question 4: Does a selected taxon differ among habitats?"), uiOutput("hypothesis_banner_taxon"),
    fluidRow(column(3, div(class = "control-card",
      selectInput("compare_rank", "Taxonomic rank", rank_choices),
      selectInput("compare_taxon", "Taxon", choices = NULL),
      checkboxGroupInput("taxon_habitats", "Habitats", habitat_levels, habitat_levels),
      checkboxGroupInput("taxon_sexes", "Include sex", sex_choices, sex_choices),
      radioButtons("taxon_test_choice", "Global test", c("Kruskal–Wallis" = "KW", "One-way ANOVA" = "ANOVA")),
      checkboxInput("taxon_pairwise", "Include pairwise comparisons", TRUE),
      actionButton("run_taxon", "Run taxon comparison", class = "btn-primary"))),
      column(9, div(class = "explanation-card",
        "This exploratory test compares relative abundance among habitats. Kruskal–Wallis uses ranks and is less sensitive to non-normal data; ANOVA compares means and assumes approximately normal, similarly variable residuals."),
        conditionalPanel("input.run_taxon == 0", div(class = "waiting-card", "Choose a taxon, groups, and test, then run the comparison.")),
        uiOutput("taxon_results_ui"))),
    div(class = "interpret-card", h3("Connect the result to your hypothesis"),
      textAreaInput("taxon_interpretation", NULL, rows = 3,
        placeholder = "Report the global result, pairwise comparisons, prevalence, and an ecological interpretation."))
  ))
)

server <- function(input, output, session) {
  all_analysis_choices <- c(
    "Alpha-diversity index + group test" = "alpha",
    "Community distance + PCoA + PERMANOVA" = "beta",
    "Relative-abundance composition graph" = "composition",
    "Selected taxon + Kruskal–Wallis/ANOVA" = "taxon")
  observeEvent(input$lab_question, {
    updateSelectInput(session, "planned_analysis", choices = all_analysis_choices, selected = character(0))
  }, ignoreInit = FALSE)
  output$analysis_feedback <- renderUI({
    req(input$planned_analysis)
    if (identical(input$planned_analysis, input$lab_question))
      tags$span(strong("Good match. "), "This analysis addresses your selected ecological question.")
    else tags$span(strong("Try again. "), "The response measured by this analysis does not match the selected question.")
  })
  observeEvent(input$begin_analysis, {
    if (!identical(input$planned_analysis, input$lab_question)) {
      showNotification("Choose an analysis that matches the ecological question.", type = "warning")
      return()
    }
    if (!nzchar(trimws(input$hypothesis))) {
      showNotification("Write a hypothesis before beginning the analysis.", type = "warning")
      return()
    }
    destination <- c(alpha = "Alpha diversity", beta = "Beta diversity", composition = "Composition", taxon = "Compare a taxon")[[input$lab_question]]
    updateNavbarPage(session, "main_nav", selected = destination)
  })

  hypothesis_ui <- function() {
    text <- trimws(if (is.null(input$hypothesis)) "" else input$hypothesis)
    if (!nzchar(text)) div(class = "hypothesis-banner", strong("Hypothesis not yet entered. "), "Return to Start and write one before interpreting results.")
    else div(class = "hypothesis-banner", strong("Your hypothesis: "), text)
  }
  output$hypothesis_banner <- renderUI(hypothesis_ui())
  output$hypothesis_banner_beta <- renderUI(hypothesis_ui())
  output$hypothesis_banner_composition <- renderUI(hypothesis_ui())
  output$hypothesis_banner_taxon <- renderUI(hypothesis_ui())

  alpha_all <- reactive({
    d <- diversity_values(counts)
    data.frame(SAMPLE_ID = colnames(counts), metadata, d, check.names = FALSE)
  })

  output$depth_plot <- renderPlot({
    req(input$data_habitats, input$data_sexes); d <- alpha_all()
    d <- d[d$HABITAT %in% input$data_habitats & d$SEX %in% input$data_sexes, ]
    d$READS <- if ("ORIGINAL_LIBRARY_SIZE" %in% names(d)) d$ORIGINAL_LIBRARY_SIZE else colSums(counts)[d$SAMPLE_ID]
    ggplot(d, aes(HABITAT, READS, fill = HABITAT)) + geom_boxplot(outlier.shape = NA, alpha = 0.8) +
      geom_jitter(width = 0.12, alpha = 0.5, size = 1.5) + scale_fill_manual(values = habitat_colors, drop = FALSE) +
      scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
      labs(x = NULL, y = "Original library size (log scale)") + theme_classic(base_size = 14) + theme(legend.position = "none")
  })
  output$sample_summary <- renderTable({
    req(input$data_habitats, input$data_sexes); d <- alpha_all()
    d <- d[d$HABITAT %in% input$data_habitats & d$SEX %in% input$data_sexes, ]
    d$ORIGINAL <- if ("ORIGINAL_LIBRARY_SIZE" %in% names(d)) d$ORIGINAL_LIBRARY_SIZE else colSums(counts)[d$SAMPLE_ID]
    d$RETAINED <- if ("RETAINED_LIBRARY_SIZE" %in% names(d)) d$RETAINED_LIBRARY_SIZE else colSums(counts)[d$SAMPLE_ID]
    do.call(rbind, lapply(split(d, droplevels(d$HABITAT)), function(z) data.frame(Habitat = as.character(z$HABITAT[1]),
      Samples = nrow(z), `Median original reads` = round(median(z$ORIGINAL)),
      `Median retained reads` = round(median(z$RETAINED)), `Median genera` = round(median(z$Richness), 1), check.names = FALSE)))
  }, striped = TRUE)

  output$alpha_metric_explanation <- renderUI(metric_text[[input$alpha_metric]])
  alpha_result <- eventReactive(input$run_alpha, {
    mode <- input$alpha_group; sexes <- if (mode == "Habitat") input$alpha_sexes else sex_choices
    habitats <- if (mode == "Habitat") input$alpha_habitats else input$alpha_one_habitat
    index <- analysis_index(mode, habitats, sexes, input$alpha_one_habitat)
    d <- alpha_all()[index, ]; d$GROUP <- if (mode == "Habitat") droplevels(d$HABITAT) else droplevels(d$SEX)
    validate(need(nlevels(d$GROUP) >= 2, "Choose at least two groups with samples.")); value <- d[[input$alpha_metric]]
    global <- if (input$alpha_test_choice == "KW") kruskal.test(value ~ d$GROUP) else summary(aov(value ~ d$GROUP))[[1]]
    pairwise <- NULL
    if (input$alpha_pairwise) pairwise <- if (input$alpha_test_choice == "KW")
      pairwise.wilcox.test(value, d$GROUP, p.adjust.method = "holm", exact = FALSE) else
      pairwise.t.test(value, d$GROUP, p.adjust.method = "holm", pool.sd = TRUE)
    list(data = d, value = value, global = global, pairwise = pairwise)
  })
  output$alpha_results_ui <- renderUI({ req(input$run_alpha > 0); tagList(
    div(class = "plot-card", plotOutput("alpha_plot", height = 440)),
    fluidRow(column(6, div(class = "table-card", h3("Summary"), tableOutput("alpha_summary"),
      p(class = "small-note", "N = sample count; median = middle value; mean = arithmetic average; IQR = spread of the middle 50%."))),
      column(6, div(class = "result-card", h3("Global test"), verbatimTextOutput("alpha_test")))),
    conditionalPanel("input.alpha_pairwise", div(class = "table-card", h3("Pairwise comparisons (Holm-adjusted)"), tableOutput("alpha_pairwise_table")))
  )})
  output$alpha_plot <- renderPlot({
    r <- alpha_result(); r$data$VALUE <- r$value
    ggplot(r$data, aes(GROUP, VALUE, fill = GROUP)) + geom_boxplot(outlier.shape = NA, alpha = 0.82) +
      geom_jitter(width = 0.12, alpha = 0.55, size = 1.7) + labs(x = NULL, y = input$alpha_metric) +
      theme_classic(base_size = 14) + theme(legend.position = "none")
  })
  output$alpha_summary <- renderTable({
    r <- alpha_result(); do.call(rbind, lapply(split(r$value, r$data$GROUP), function(z)
      data.frame(N = length(z), Median = median(z), Mean = mean(z), IQR = IQR(z))))
  }, rownames = TRUE, digits = 3, striped = TRUE)
  output$alpha_test <- renderPrint(print(alpha_result()$global))
  output$alpha_pairwise_table <- renderTable({ req(alpha_result()$pairwise); pairwise_matrix_table(alpha_result()$pairwise) }, digits = 4, striped = TRUE)

  output$distance_explanation <- renderUI({
    if (input$distance_choice == "Bray-Curtis") tags$span(strong("Bray–Curtis"), " compares relative abundances. Shared dominant taxa make samples more similar; abundance differences increase distance.")
    else tags$span(strong("Jaccard"), " uses presence and absence only. Every detected genus has equal weight, regardless of abundance.")
  })
  beta_result <- eventReactive(input$run_beta, {
    mode <- input$beta_group; sexes <- if (mode == "Habitat") input$beta_sexes else sex_choices
    habitats <- if (mode == "Habitat") input$beta_habitats else input$beta_one_habitat
    index <- analysis_index(mode, habitats, sexes, input$beta_one_habitat)
    group <- if (mode == "Habitat") droplevels(metadata$HABITAT[index]) else droplevels(metadata$SEX[index])
    validate(need(nlevels(group) >= 2, "Choose at least two groups with samples."))
    d <- distance_matrix(counts[, index, drop = FALSE], input$distance_choice)
    fit <- cmdscale(d, k = 2, eig = TRUE, add = TRUE); positive <- fit$eig[fit$eig > 0]
    points <- data.frame(PCoA1 = fit$points[, 1], PCoA2 = fit$points[, 2], GROUP = group)
    withProgress(message = "Running fixed permutation tests", value = 0.2, {
      perma <- permanova_stats(d, group, 999); incProgress(0.35)
      dispersion <- dispersion_stats(d, group, 999); incProgress(0.35)
      pairwise <- if (input$beta_pairwise && nlevels(group) > 2) pairwise_beta_tests(d, group, 199) else NULL
      incProgress(0.1)
    })
    list(points = points, axis = round(100 * fit$eig[1:2] / sum(positive), 1),
      permanova = perma, dispersion = dispersion, pairwise = pairwise)
  })
  output$beta_results_ui <- renderUI({ req(input$run_beta > 0); tagList(
    div(class = "plot-card", plotOutput("pcoa_plot", height = 480),
      div(class = "figure-note", "Each point is a community. Nearby points are more similar. Axis percentages show how much corrected distance variation each axis represents.")),
    fluidRow(column(6, div(class = "table-card", h3("PERMANOVA: group centroids"), tableOutput("permanova_table"))),
      column(6, div(class = "table-card", h3("Dispersion: within-group variability"), tableOutput("dispersion_table")))),
    conditionalPanel("input.beta_pairwise", div(class = "table-card", h3("Pairwise tests"),
      p(class = "small-note", "Pairwise tests use 199 permutations; p-values are Holm-adjusted."), tableOutput("beta_pairwise_table"))),
    div(class = "callout", strong("Interpret together: "), "PERMANOVA tests group centroids. Significant dispersion means within-group variability also differs.")
  )})
  output$pcoa_plot <- renderPlot({
    r <- beta_result(); ggplot(r$points, aes(PCoA1, PCoA2, color = GROUP)) + geom_point(size = 2.7, alpha = 0.78) +
      labs(x = paste0("PCoA axis 1 (", r$axis[1], "%)"), y = paste0("PCoA axis 2 (", r$axis[2], "%)"), color = "Group") + theme_classic(base_size = 14)
  })
  output$permanova_table <- renderTable(beta_result()$permanova, digits = 4, striped = TRUE)
  output$dispersion_table <- renderTable(beta_result()$dispersion, digits = 4, striped = TRUE)
  output$beta_pairwise_table <- renderTable({
    x <- beta_result()$pairwise; req(x); x$PERMANOVA_p <- p.adjust(x$PERMANOVA_p, "holm")
    x$Dispersion_p <- p.adjust(x$Dispersion_p, "holm"); names(x)[3:4] <- c("PERMANOVA adjusted p", "Dispersion adjusted p"); x
  }, digits = 4, striped = TRUE)

  composition_result <- eventReactive(input$run_composition, {
    req(input$composition_habitats, input$composition_sexes)
    index <- metadata$HABITAT %in% input$composition_habitats & metadata$SEX %in% input$composition_sexes
    validate(need(any(index), "No samples match these habitat and sex filters."))
    rank_counts <- rank_tables[[input$tax_rank]][, index, drop = FALSE]; sub_meta <- metadata[index, , drop = FALSE]
    if (input$composition_view == "Mean") composition_summary(rank_counts, sub_meta, input$top_n)
    else sample_composition(rank_counts, sub_meta, input$top_n)
  })
  output$composition_results_ui <- renderUI({ req(input$run_composition > 0); tagList(
    div(class = "plot-card", plotOutput("composition_plot", height = 520)),
    conditionalPanel("input.composition_view == 'Mean'", div(class = "table-card", h3("Mean relative abundance (%)"), tableOutput("composition_table")))
  )})
  output$composition_plot <- renderPlot({
    d <- composition_result()
    if (input$composition_view == "Mean") ggplot(d, aes(HABITAT, ABUNDANCE, fill = TAXON)) +
      geom_col(color = "white", linewidth = 0.2) + scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1), expand = c(0, 0)) +
      labs(x = NULL, y = "Mean relative abundance", fill = input$tax_rank) + theme_classic(base_size = 14)
    else ggplot(d, aes(SAMPLE, ABUNDANCE, fill = TAXON)) + geom_col(width = 1) +
      facet_grid(~ HABITAT, scales = "free_x", space = "free_x") +
      scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1), expand = c(0, 0)) +
      labs(x = "Individual samples", y = "Relative abundance", fill = input$tax_rank) + theme_classic(base_size = 13) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = grid::unit(0.15, "lines"))
  })
  output$composition_table <- renderTable({
    d <- composition_result(); wide <- reshape(transform(d, Percent = round(100 * ABUNDANCE, 1))[, c("HABITAT", "TAXON", "Percent")],
      idvar = "TAXON", timevar = "HABITAT", direction = "wide")
    names(wide) <- sub("Percent\\.", "", names(wide)); wide
  }, digits = 1, striped = TRUE)

  observeEvent(input$compare_rank, {
    taxa <- rownames(rank_tables[[input$compare_rank]])
    updateSelectInput(session, "compare_taxon", choices = sort(taxa), selected = taxa[1])
  }, ignoreInit = FALSE)
  taxon_result <- eventReactive(input$run_taxon, {
    req(input$compare_taxon, input$taxon_habitats, input$taxon_sexes)
    rank_counts <- rank_tables[[input$compare_rank]]
    index <- metadata$HABITAT %in% input$taxon_habitats & metadata$SEX %in% input$taxon_sexes
    validate(need(any(index), "No samples match these habitat and sex filters."))
    sub_meta <- metadata[index, , drop = FALSE]
    abundance <- rank_counts[input$compare_taxon, index] / colSums(rank_counts[, index, drop = FALSE])
    group <- droplevels(sub_meta$HABITAT); validate(need(nlevels(group) >= 2, "Choose at least two habitats with samples."))
    global <- if (input$taxon_test_choice == "KW") kruskal.test(abundance ~ group) else summary(aov(abundance ~ group))[[1]]
    pairwise <- NULL
    if (input$taxon_pairwise) pairwise <- if (input$taxon_test_choice == "KW")
      pairwise.wilcox.test(abundance, group, p.adjust.method = "holm", exact = FALSE) else
      pairwise.t.test(abundance, group, p.adjust.method = "holm", pool.sd = TRUE)
    plot_data <- data.frame(HABITAT = group, ABUNDANCE = abundance)
    prevalence <- do.call(rbind, lapply(split(abundance, group), function(z) data.frame(N = length(z),
      `Prevalence (%)` = 100 * mean(z > 0), `Median abundance (%)` = 100 * median(z),
      `Mean abundance (%)` = 100 * mean(z), check.names = FALSE)))
    list(data = plot_data, global = global, pairwise = pairwise, prevalence = prevalence)
  })
  output$taxon_results_ui <- renderUI({ req(input$run_taxon > 0); tagList(
    div(class = "plot-card", plotOutput("taxon_plot", height = 430)),
    fluidRow(column(6, div(class = "table-card", h3("Prevalence and abundance"), tableOutput("taxon_prevalence"))),
      column(6, div(class = "result-card", h3("Global test"), verbatimTextOutput("taxon_test")))),
    conditionalPanel("input.taxon_pairwise", div(class = "table-card", h3("Pairwise comparisons (Holm-adjusted)"), tableOutput("taxon_pairwise_table"))),
    div(class = "callout", strong("Statistical caution: "),
      "relative abundances are bounded, often contain zeros, and are compositional. Treat this classroom test as exploratory."))
  })
  output$taxon_plot <- renderPlot({
    r <- taxon_result(); ggplot(r$data, aes(HABITAT, 100 * ABUNDANCE, fill = HABITAT)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.82) + geom_jitter(width = 0.12, alpha = 0.55, size = 1.6) +
      scale_fill_manual(values = habitat_colors, drop = FALSE) +
      labs(x = NULL, y = paste(input$compare_taxon, "relative abundance (%)")) + theme_classic(base_size = 14) + theme(legend.position = "none")
  })
  output$taxon_prevalence <- renderTable(taxon_result()$prevalence, rownames = TRUE, digits = 2, striped = TRUE)
  output$taxon_test <- renderPrint(print(taxon_result()$global))
  output$taxon_pairwise_table <- renderTable({ req(taxon_result()$pairwise); pairwise_matrix_table(taxon_result()$pairwise) }, digits = 4, striped = TRUE)
}

shinyApp(ui, server)
