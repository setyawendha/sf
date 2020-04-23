# see https://docs.google.com/presentation/d/1Hl4KapfAENAOf4gv-pSngKwvS_jwNVHRPZTTDzXXn6Q/view?pli=1#slide=id.i0
# and the S2 package:
# https://cran.r-project.org/package=s2
# https://github.com/spatstat/s2

s2_p4s = "+proj=geocent +a=6371000 +b=6371000"
earth_radius = 6371000.0

#nocov start

#' functions for spherical geometry, using s2 package
#' 
#' functions for spherical geometry, using the s2 package based on the google s2 library
#' @name s2 
#' @param to_p4s target CRS; defaults to "+proj=geocent +a=6371000 +b=6371000"
#' @export
#' @details \code{st_as_s2} converts an \code{sf} POLYGON object into a form readable by \code{s2}.
st_as_s2 = function(x, to_p4s = s2_p4s) {
	# geocentric, spherical:
	geom = st_transform(st_geometry(x), st_crs(to_p4s))
	g = opp_sfc(geom, as.numeric(0.0), 0L, NA_crs_, 1L) # computes unit-length XYZ vectors
	if (inherits(geom, "sfc_MULTIPOLYGON")) # unlist: S2 sorts out what are holes
		g = lapply(g, unlist, recursive = FALSE)
	else if (!inherits(geom, "sfc_POLYGON"))
		stop("only objects of class sfc_MULTIPOLYGON or sfc_POLYGON accepted")
	g
}

load_s2 = function() {
	if (! requireNamespace("s2", quietly = TRUE))
		stop("package s2 required, please install it first")
}

#' @name s2
#' @export
#' @param x object of class \code{sf}, or \code{sfc}
#' @param ... passed on to \link{st_transform}
#' @param crs object of class \code{crs}
#' @param radius numeric; Earth radius in m; defaults to 6371000.0
st_as_sfc.S2Polygon = function(x, ..., crs = st_crs(4326), radius = earth_radius) {
	# close all loops:
	loops = lapply(x$loops, function(L) rbind(L, L[1,]) * radius)
	loops[x$holes] = lapply(loops[x$holes], function(L) L[nrow(L):1, ])
	w = which(! x$holes)
	splt = rep(seq_along(w), diff(c(w, length(x$holes) + 1)))
	p = if (length(splt) > 1)
		st_multipolygon(split(loops, splt))
	else
		st_polygon(loops)
	st_zm(st_transform(st_sfc(p, crs = st_crs(s2_p4s)), crs, ...))
}

#' @export
st_as_sfc.S2Polygons = function(x, ..., crs = st_crs(4326), radius = earth_radius) {
	structure(do.call(c, lapply(x, st_as_sfc, crs = crs, radius = radius)), n_empty = attr(x, "n_empty"))
}

#' @export
st_as_sfc.S2Points = function(x, ..., crs = st_crs(4326), radius = earth_radius) {
	x = x * radius
	st_zm(st_transform(
		st_sfc(lapply(1:nrow(x), function(i) st_point(x[i,])), crs = st_crs(s2_p4s)),
		crs))
}

#' @name s2
#' @export
#' @param y object of class \code{sf}, or \code{sfc}
s2_intersection = function(x, y) {
	load_s2()
	stopifnot(st_crs(x) == st_crs(y))
	st_as_sfc(s2::S2Polygons_intersection(st_as_s2(x), st_as_s2(y)), crs = st_crs(x))
}

#' @name s2
#' @export
s2_intersects = function(x, y = x) {
	load_s2()
	stopifnot(st_crs(x) == st_crs(y))
	s2::S2Polygons_intersect(st_as_s2(x), st_as_s2(y))
}

#' @name s2
#' @export
#' @details \code{s2_centroid} computes the spherical centroid of a set of polygons
#' @examples
#' demo(nc, echo = FALSE, ask = FALSE)
#' s2_centroid(nc)
#' pts = rbind(c(0,80), c(120,80), c(240,80), c(0,80))
#' pole = st_sfc(st_polygon(list(pts)), crs = 4326)
#' s2_centroid(pole)
s2_centroid = function(x) {
	load_s2()
	st_as_sfc(structure(s2::S2Polygons_centroid(st_as_s2(x)), class = "S2Points"),
		crs = st_crs(x))
}

#' @name s2
#' @export
#' @details \code{s2_area} computes the area of a set of polygons, as a fraction of 4 * pi.
#' @examples
#' demo(nc, echo = FALSE, ask = FALSE)
#' s2_area(nc)
s2_area = function(x) {
	load_s2()
	s2::S2Polygons_area(st_as_s2(x)) * earth_radius^2
}

#nocov end
