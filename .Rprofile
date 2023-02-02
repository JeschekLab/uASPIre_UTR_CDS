# set repositories to fixed server
options(repos = structure(c(CRAN = 'https://stat.ethz.ch/CRAN/')))

# set theme for plotting
theme_SH <- function(font_size = 14, font_small = 12) {
	theme_bw() +
	theme(
		# text
		text = element_text(
			family = 'ArialMT',
      face = 'plain',
      colour = '#000000',
      size = font_size,
      hjust = 0.5,
      vjust = 0.5,
      angle = 0,
      lineheight = 0.9,
      margin = margin(),
      debug = FALSE),
		# plot area
		plot.background = element_blank(),
		# panel
		panel.background = element_blank(),
		panel.border = element_rect(
			fill = NA,
			colour = '#000000',
			linewidth = 0.46875),
		panel.ontop = FALSE,
		panel.grid.major = element_line(
			linewidth = 0.234375),
		panel.grid.minor = element_blank(),
		# axis text, title and size
		axis.text = element_text(
			color = '#000000',
			family = 'ArialMT',
			size = font_small),
		axis.title = element_text(
			colour = '#000000',
			family = 'ArialMT',
			size = font_size),
		axis.ticks = element_line(
			colour = '#000000',
			linewidth = 0.46875),
		axis.ticks.length = unit(2.51, 'points'),
		# legend
		legend.background = element_blank(),
    legend.spacing = unit(font_size, 'pt'),
    legend.spacing.x = NULL,
    legend.spacing.y = NULL,
    legend.margin = margin(0, 0, 0, 0),
    legend.key = element_blank(),
    legend.key.size = unit(1.1 * font_size, 'pt'),
    legend.key.height = NULL,
    legend.key.width = NULL,
    legend.text = element_text(size = rel(font_small/font_size)),
    legend.text.align = NULL,
    legend.title = element_text(hjust = 0),
    legend.title.align = NULL,
    legend.position = 'right',
    legend.direction = NULL,
    legend.justification = c('left', 'center'),
    legend.box = NULL,
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.background = element_blank(),
    legend.box.spacing = unit(font_size, 'pt')
	)
}

# show function
HEAD <- function(x) {as.data.frame(head(x))}
TAIL <- function(x) {as.data.frame(tail(x))}
