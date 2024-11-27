#include <gimp-2.0/libgimp/gimp.h>
#include <gimp-2.0/libgimp/gimpui.h>
#include <math.h>
#include <stdlib.h>

guint8 get_pixel_intensity(gboolean is_rgb, const guchar* pixel) {
	if(is_rgb)
		return ROUND(0.2126 * pixel[0] + 0.7154 * pixel[1] + 0.0722 * pixel[2]);
	return *pixel;
}

guint8 apply_gamma_adjusting(guint8 intensity) {
	gdouble scaled = (gdouble)intensity / 255.0;
	if(scaled <= 0.0031308)
		return ROUND(12.92 * scaled * 255.0);
	return ROUND(255.0 * (1.055 * pow(scaled, 1.0 / 2.4) - 0.055));	
}

void get_grayscaled_row(gboolean is_rgb, gint bpp, gint width, guchar* inrow, guchar* outrow) {
	for (guint i = 0; i < width * bpp; i += bpp) {
		guint8 intensity = get_pixel_intensity(is_rgb, inrow + i);
		intensity = apply_gamma_adjusting(intensity);
		if(is_rgb) {
			outrow[i] = intensity;
			outrow[i + 1] = intensity;
			outrow[i + 2] = intensity;
		} else 
			outrow[i] = intensity;
	}
}

void convert_grayscale(GimpDrawable* drawable) {
	guint width = drawable->width;
	guint height = drawable->height;
	gint bpp = drawable->bpp;
	gint32 id = drawable->drawable_id;
		
	GimpPixelRgn rgn_in, rgn_out;
	guchar *inrow, *outrow;

	gimp_pixel_rgn_init(&rgn_in, drawable, 0, 0, width, height, FALSE, FALSE);
	gimp_pixel_rgn_init(&rgn_out, drawable, 0, 0, width, height, TRUE, TRUE);

	inrow = g_new(guchar, bpp * width);
	outrow = g_new(guchar, bpp * width);

	for(gint y = 0; y < height; ++y) {
		gimp_pixel_rgn_get_row(&rgn_in, inrow, 0, y, width);
		memcpy(outrow, inrow, bpp * width);
		get_grayscaled_row(gimp_drawable_is_rgb(id), bpp, width, inrow, outrow);
		gimp_pixel_rgn_set_row(&rgn_out, outrow, 0, y, width);
	}

	g_free(inrow);
	g_free(outrow);
	gimp_drawable_flush(drawable);
	gimp_drawable_merge_shadow(id, TRUE);
	gimp_drawable_update(id, 0, 0, width, height);
}

typedef struct {
	gdouble* data;
	gint w, h;
	gint bpp;
	gint rowstride;
	gboolean has_alpha;
} AdjustedPixelRgn;


AdjustedPixelRgn* init_rgn(gint w, gint h, gint bpp, gboolean has_alpha) {
	AdjustedPixelRgn* rgn = g_new(AdjustedPixelRgn, 1);
	rgn->w = w;
	rgn->h = h;
	rgn->bpp = bpp;
	rgn->rowstride = w * bpp;
	rgn->has_alpha = has_alpha;
	rgn->data = g_new(gdouble, w * h * bpp);
	return rgn;
}

AdjustedPixelRgn* copy(const AdjustedPixelRgn* src) {
	AdjustedPixelRgn* copy = init_rgn(src->w, src->h, src->bpp, src->has_alpha);
	memcpy(copy->data, src->data, src->w * src->h * src->bpp * sizeof(gdouble));
	return copy;
}

void get_pixel_rgn_row(const AdjustedPixelRgn* rgn, gdouble* buf, gint x, gint y, gint w) {
	memcpy(buf, rgn->data + y * rgn->rowstride + x * rgn->bpp, w * rgn->bpp * sizeof(gdouble));
}

void get_pixel_rgn_rect(const AdjustedPixelRgn* rgn, gdouble* buf, gint x, gint y, gint w, gint h) {
	for(gint row = 0; row < h; ++row) {
		get_pixel_rgn_row(rgn, buf + row * w * rgn->bpp, x, y + row, w);
	}
}

void get_pixel_rgn_pixel(const AdjustedPixelRgn* rgn, gdouble* buf, gint x, gint y) {
	get_pixel_rgn_row(rgn, buf, x, y, 1);
}

void get_pixel_rgn_col(const AdjustedPixelRgn* rgn, gdouble* buf, gint x, gint y, gint h) {
	get_pixel_rgn_rect(rgn, buf, x, y, 1, h);
}

void set_pixel_rgn_row(AdjustedPixelRgn* rgn, const gdouble* buf, gint x, gint y, gint w) {
	memcpy(rgn->data + y * rgn->rowstride + x * rgn->bpp, buf, w * rgn->bpp * sizeof(gdouble));
}

void set_pixel_rgn_pixel(AdjustedPixelRgn* rgn, const gdouble* buf, gint x, gint y) {
	set_pixel_rgn_row(rgn, buf, x, y, 1);	
}

void set_pixel_rgn_rect(AdjustedPixelRgn* rgn, const gdouble* buf, gint x, gint y, gint w, gint h) {
	for(gint row = 0; row < h; ++row) {
		set_pixel_rgn_row(rgn, buf + row * w * rgn->bpp, x, y + row, w);
	}
}

void set_pixel_rgn_col(AdjustedPixelRgn* rgn, const gdouble* buf, gint x, gint y, gint h) {
	set_pixel_rgn_rect(rgn, buf, x, y, 1, h);
}

void free_rgn(AdjustedPixelRgn* rgn) {
	g_free(rgn->data);
	g_free(rgn);
}

void convert_pixel_to_double(const guchar* src, gdouble* dest, guint len) {
	for(guint i = 0; i < len; ++i) 
		dest[i] = (gdouble)src[i];
}

void adjust_top_corners(AdjustedPixelRgn* to, GimpPixelRgn* from, guint8 r) {
	gint w = from->w, h = from->h, bpp = from->bpp;
	guchar* top_left = g_new(guchar, bpp);
	guchar* top_right = g_new(guchar, bpp);
	gimp_pixel_rgn_get_pixel(from, top_left, 0, 0);
	gimp_pixel_rgn_get_pixel(from, top_right, w - 1, 0);
	gdouble* converted_top_left = g_new(gdouble, bpp);
	gdouble* converted_top_right = g_new(gdouble, bpp);
	convert_pixel_to_double(top_left, converted_top_left, bpp);
	convert_pixel_to_double(top_right, converted_top_right, bpp);
	
	for(gint y = 0; y < r; ++y) {
		for(gint x = 0; x < r; ++x) {
			set_pixel_rgn_pixel(to, converted_top_left, x, y);
			set_pixel_rgn_pixel(to, converted_top_right, w + r + x, y);
		}
	}
	
	g_free(top_right);
	g_free(top_left);
	g_free(converted_top_right);
	g_free(converted_top_left);
}

void adjust_bottom_corners(AdjustedPixelRgn* to, GimpPixelRgn* from, guint8 r) {
	gint w = from->w, h = from->h, bpp = from->bpp;
	guchar* bottom_left = g_new(guchar, bpp);
	guchar* bottom_right = g_new(guchar, bpp);
	gimp_pixel_rgn_get_pixel(from, bottom_left, 0, h - 1);
	gimp_pixel_rgn_get_pixel(from, bottom_right, w - 1, h - 1);
	gdouble* converted_bottom_left = g_new(gdouble, bpp);
	gdouble* converted_bottom_right = g_new(gdouble, bpp);
	convert_pixel_to_double(bottom_left, converted_bottom_left, bpp);
	convert_pixel_to_double(bottom_right, converted_bottom_right, bpp);

	for(gint y = 0; y < r; ++y) {
		for(gint x = 0; x < r; ++x) {
			set_pixel_rgn_pixel(to, converted_bottom_left, x, h + r + y);
			set_pixel_rgn_pixel(to, converted_bottom_right, w + r + x, h + r + y);
		}
	}
		
	g_free(bottom_right);
	g_free(bottom_left);
	g_free(converted_bottom_right);
	g_free(converted_bottom_left);
}

void adjust_top(AdjustedPixelRgn* to, GimpPixelRgn* from, guint8 r) {
	gint w = from->w, h = from->h, bpp = from->bpp;
	guchar* row = g_new(guchar, w * bpp);
	gimp_pixel_rgn_get_row(from, row, 0, 0, w);
	gdouble* converted_row = g_new(gdouble, w * bpp);
	convert_pixel_to_double(row, converted_row, w * bpp);
	
	for(gint y = 0; y < r; ++y) {
		set_pixel_rgn_row(to, converted_row, r, y, w);
	}
	g_free(row);
	g_free(converted_row);

	adjust_top_corners(to, from, r);
}

void adjust_bottom(AdjustedPixelRgn* to, GimpPixelRgn* from, guint8 r) {
	gint w = from->w, h = from->h, bpp = from->bpp;
	guchar* row = g_new(guchar, w * bpp);
	gimp_pixel_rgn_get_row(from, row, w - 1, h - 1, w);
	gdouble* converted_row = g_new(gdouble, w * bpp);
	convert_pixel_to_double(row, converted_row, w * bpp);
	for(gint y = 0; y < r; ++y) {
		set_pixel_rgn_row(to, converted_row, r, r + h + y, w);
	}
	g_free(row);
	g_free(converted_row);

	adjust_bottom_corners(to, from, r);
}

void adjust_sides(AdjustedPixelRgn* to, GimpPixelRgn* from, guint8 r) {
	gint w = from->w, h = from->h, bpp = from->bpp;
	guchar* left_col = g_new(guchar, h * bpp);
	guchar* right_col = g_new(guchar, h * bpp);
	gimp_pixel_rgn_get_col(from, left_col, 0, 0, h);
	gimp_pixel_rgn_get_col(from, right_col, w - 1, 0, h);
	gdouble* converted_left_col = g_new(gdouble, h * bpp);
	gdouble* converted_right_col = g_new(gdouble, h * bpp);
	convert_pixel_to_double(left_col, converted_left_col, h * bpp);
	convert_pixel_to_double(right_col, converted_right_col, h * bpp);
	
	for(gint x = 0; x < r; ++x) {
		set_pixel_rgn_col(to, converted_left_col, x, r, h);
		set_pixel_rgn_col(to, converted_right_col, w + r + x, r, h);
	}
	g_free(left_col);
	g_free(right_col);
	g_free(converted_left_col);
	g_free(converted_right_col);
}

void adjust_inner(AdjustedPixelRgn* to, GimpPixelRgn* from, guint8 r) {
	gint w = from->w, h = from->h, bpp = from->bpp;
	
	guchar* buf = g_new(guchar, w * h * bpp);
	gimp_pixel_rgn_get_rect(from, buf, 0, 0, w, h);
	gdouble* converted_buf = g_new(gdouble, w * h * bpp);
	convert_pixel_to_double(buf, converted_buf, w * h * bpp);
	set_pixel_rgn_rect(to, converted_buf, r, r, w, h);
	g_free(buf);
	g_free(converted_buf);
}

AdjustedPixelRgn* adjust_pixel_rgn(GimpPixelRgn* in, guint8 kernel_size, gboolean has_alpha) {
	guint8 r = kernel_size / 2;
	AdjustedPixelRgn* rgn = init_rgn(in->w + 2 * r, in->h + 2 * r, in->bpp, has_alpha);
	adjust_top(rgn, in, r);
	adjust_bottom(rgn, in, r);
	adjust_sides(rgn, in, r);
	adjust_inner(rgn, in, r);
	return rgn;		
}

guchar fit(gdouble intensity) {
	if (intensity < 0.0)
		return (guchar)0;
	if (intensity > 255.0)
		return (guchar)255;
	return (guchar)ROUND(intensity);		
}

void fit_buffer(const gdouble* buf, guchar* out, guint len) {
	for(guint i = 0; i < len; ++i) 
		out[i] = fit(buf[i]);
}

void inverse_adjust_pixel_rgn(AdjustedPixelRgn* rgn, GimpPixelRgn* out, guint8 kernel_size) {
	gint w = out->w, h = out->h, bpp = out->bpp;
	guint8 r = kernel_size / 2;
	
	guchar* out_buf = g_new(guchar, w * h * bpp);
	gdouble* in_buf = g_new(gdouble, w * h * bpp);
	get_pixel_rgn_rect(rgn, in_buf, r, r, w, h);
	fit_buffer(in_buf, out_buf, w * h * bpp);
	gimp_pixel_rgn_set_rect(out, out_buf, 0, 0, w, h);
	g_free(in_buf);
	g_free(out_buf);
}

void get_convolution_window_for_pixel(const AdjustedPixelRgn* rgn, guint8 kernel_size, gint x, gint y, gdouble* window) {
	guint8 r = kernel_size / 2;
	get_pixel_rgn_rect(rgn, window, x - r, y - r, kernel_size, kernel_size);
}

void calculate_weighted_sum(const gdouble* values, const gdouble* weights, gdouble* res, const gint len, const gint bpp, gboolean has_alpha) {
	gint max_channel = bpp;
	if(has_alpha) 
		res[--max_channel] = values[max_channel];
	gdouble sum = 0.0;
	for(gint i = 0; i < len; ++i) 
		sum += values[i * bpp] * weights[i];
	for(gint channel = 0; channel < max_channel; ++channel)
		res[channel] = sum;
}

AdjustedPixelRgn* get_convolution_result(const AdjustedPixelRgn* in, const gdouble* kernel, guint8 kernel_size) {
	AdjustedPixelRgn* result = copy(in);
	gdouble* window = g_new(gdouble, kernel_size * kernel_size * in->bpp);
	gdouble* pixel = g_new(gdouble, in->bpp);
	guint8 r = kernel_size / 2;
	gint w = in->w - 2 * r, bpp = in->bpp;
	gdouble* row = g_new(gdouble, w * bpp);
	
	for(gint y = r; y < in->h - r; ++y) {
		gdouble* row_ptr = row;
		for(gint x = r; x < in->w - r; ++x, row_ptr += bpp) {
			get_convolution_window_for_pixel(in, kernel_size, x, y, window);
			calculate_weighted_sum(window, kernel, row_ptr, kernel_size * kernel_size, in->bpp, in->has_alpha);
		}
		set_pixel_rgn_row(result, row, r, y, w);
	}
	
	g_free(row);
	g_free(pixel);
	g_free(window);	
	return result;
}

guint8 calculate_gaussian_kernel_size(gdouble sigma) {
	return 2 * ((guint8)ceil(3 * sigma)) + 1;
}

void calculate_gaussian_kernel(gdouble sigma, guint8 size, gdouble* kernel) {
	guint8 r = size / 2;
	gdouble sum = 0.0;
	for(gint y = 0; y < size; ++y) {
		for(gint x = 0; x < size; ++x) {
			gdouble pos = (x - r) * (x - r) + (y - r) * (y - r);
			gdouble power = - pos / (2 * sigma * sigma);
			sum += exp(power);
			*(kernel + y * size + x) = exp(power);
		}
	}
	
	for(gint i = 0; i < size * size; ++i) {
		kernel[i] /= sum;
	}
}

typedef struct {
	AdjustedPixelRgn* magnitude;
	AdjustedPixelRgn* angle;
} ImageGradient;

gdouble get_fixed_direction(gdouble angle) {
	if (angle < 22.5 && angle >=0 || angle <= 180 && angle >= 157.5)
		return 0.0;
	if(angle >= 22.5 && angle < 67.5)
		return 45.0;
	if(angle >= 67.5 && angle < 112.5)
		return 90.0;
	if(angle >= 112.5 && angle < 157.5)
		return 135.0;			
}

gdouble radian_to_degree_fixed(gdouble radians) {
	gdouble degrees = radians * 180.0 / M_PI;
	if(degrees < 0.0)
		degrees += 180.0;
	return get_fixed_direction(degrees);	
}

void get_pixel_gradient(const gdouble* g_x, const gdouble* g_y, gdouble* mag, gdouble* angle, gint bpp, gboolean has_alpha, gdouble* max) {
	if(has_alpha)
		mag[--bpp] = g_x[bpp];
	gdouble magnitude = hypot(*g_x, *g_y);
	if(magnitude > *max)
		*max = magnitude;
	gdouble theta = atan2(*g_y, *g_x);	
	*angle = radian_to_degree_fixed(theta);
	for(gint i = 0; i < bpp; ++i)
		mag[i] = magnitude;
}

ImageGradient get_image_gradient(const AdjustedPixelRgn* G_x, const AdjustedPixelRgn* G_y, guint w, guint h, guint radius) {
	ImageGradient gradient;
	gint bpp = G_x->bpp;
	gradient.magnitude = init_rgn(w, h, bpp, G_x->has_alpha);
	gradient.angle = init_rgn(w, h, 1, FALSE);
	gdouble* mag = gradient.magnitude->data;
	gdouble* angle = gradient.angle->data;
	gdouble* g_x = g_new(gdouble, w * h * bpp);
	gdouble* g_y = g_new(gdouble, w * h * bpp);
	get_pixel_rgn_rect(G_x, g_x, radius, radius, w, h);
	get_pixel_rgn_rect(G_y, g_y, radius, radius, w, h);
	gdouble max = 0.0;
	
	for(gint i = 0; i < w * h; ++i, ++angle, mag += bpp) 
		get_pixel_gradient(g_x + i * bpp, g_y + i * bpp, mag, angle, bpp, G_x->has_alpha, &max);

	g_free(g_x);
	g_free(g_y);	

	return gradient;
}

gboolean is_local_maxima(const gdouble* window, gdouble angle, gint bpp) {
	gint first_offset, second_offset;
	gint center = 4 * bpp;
	if(angle == 0.0) {
		first_offset = 3 * bpp;
		second_offset = 5 * bpp;
	}
	if(angle == 45.0) {
		first_offset = 2 * bpp;
		second_offset = 6 * bpp;
	}
	if(angle == 90.0) {
		first_offset = 1 * bpp;
		second_offset = 7 * bpp;
	}
	if(angle == 135.0) {
		first_offset = 0;
		second_offset = 8 * bpp;
	}
	return (window[center] >= window[first_offset]) && (window[center] >= window[second_offset]);
}

void supress(gdouble* pixel, gint bpp, gboolean has_alpha) {
	if(has_alpha)
		--bpp;
	for(gint i = 0; i < bpp; ++i)
		pixel[i] = 0.0;	
}

gdouble classify_pixel(gdouble intensity, gdouble high_threshold, gdouble low_threshold) {
	if(intensity >= high_threshold)
		return 1.0;
	if (intensity < low_threshold)
		return 0.0;
	return -1.0;		
}

AdjustedPixelRgn* map_pixels(AdjustedPixelRgn* rgn, gdouble high_threshold, gdouble low_threshold) {
	gint w = rgn->w, h = rgn->h, bpp = rgn->bpp;
	AdjustedPixelRgn* mapping = init_rgn(w, h, 1, FALSE);
	gdouble* data_ptr = rgn->data;
	gdouble* map_ptr = mapping->data;
	for(guint i = 0; i < w * h; ++i, data_ptr += bpp, ++map_ptr){
		*map_ptr = classify_pixel(*data_ptr, high_threshold, low_threshold);
	}
	return mapping;
}

gboolean has_connection_to_strong_edge(const gdouble* mapping_window) {
	for(guint i = 0; i < 9; ++i) {
		if(mapping_window[i] == 1.0)
			return TRUE;
	}
	return FALSE;
}

void apply_hysteresis_tracking(AdjustedPixelRgn* rgn, gdouble high_threshold, gdouble low_threshold) {
	AdjustedPixelRgn* mapping = map_pixels(rgn, high_threshold, low_threshold);
	gint w = rgn->w, h = rgn->h, bpp = rgn->bpp;
	gboolean has_alpha = rgn->has_alpha;
	gdouble* map_window = g_new(gdouble, 3 * 3);

	for(guint y = 1; y < h - 1; ++y) {
		gdouble* data_ptr = rgn->data + y * rgn->rowstride + 1 * bpp;
		gdouble* map_ptr = mapping->data + y * mapping->rowstride + 1;
		for(guint x = 1; x < w - 1; ++x, ++map_ptr) {
			get_convolution_window_for_pixel(mapping, 3, x, y, map_window);
			gboolean to_supress = (*map_ptr == -1.0) || (*map_ptr == 0.0 && !has_connection_to_strong_edge(map_window));
			if(to_supress)
				supress(data_ptr, bpp, has_alpha);
		}
	}
	free_rgn(mapping);
	g_free(map_window);
}

AdjustedPixelRgn* apply_non_max_suppression(ImageGradient* gradient) {
	AdjustedPixelRgn* result = copy(gradient->magnitude);
	gint bpp = result->bpp, h = result->h - 1, w = result->w - 1;
	gdouble* window = g_new(gdouble, 3 * 3 * bpp);
	
	for(gint y = 1; y < h; ++y) {
		gdouble* mag_ptr = result->data + y * result->rowstride + 1 * bpp;
		gdouble* angle_ptr = gradient->angle->data + y * gradient->angle->rowstride + 1;
		for(gint x = 1; x < w ; ++x, mag_ptr += bpp, ++angle_ptr) {
			get_convolution_window_for_pixel(gradient->magnitude, 3, x, y, window);
			if(!is_local_maxima(window, *angle_ptr, bpp)) 
				supress(mag_ptr, bpp, result->has_alpha);
		}
	}

	g_free(window);
	return result;
}

ImageGradient apply_sobel(GimpDrawable* drawable) {
	guint width = drawable->width;
	guint height = drawable->height;
	gint bpp = drawable->bpp;
	gint32 id = drawable->drawable_id;
			
	GimpPixelRgn rgn_in, rgn_out;
		
	gimp_pixel_rgn_init(&rgn_in, drawable, 0, 0, width, height, FALSE, FALSE);
	gimp_pixel_rgn_init(&rgn_out, drawable, 0, 0, width, height, TRUE, TRUE);

	guint8 kernel_size = 3;
	gdouble sobel_x[9] = {-0.25, 0.0, 0.25, -0.5, 0.0, 0.5, -0.25, 0.0, 0.25};
	gdouble sobel_y[9] = {0.25, 0.5, 0.25, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25};
	gboolean has_alpha = gimp_drawable_has_alpha(id);

	AdjustedPixelRgn* adj_in = adjust_pixel_rgn(&rgn_in, kernel_size, has_alpha);
	AdjustedPixelRgn* G_x = get_convolution_result(adj_in, sobel_x, kernel_size);
	AdjustedPixelRgn* G_y = get_convolution_result(adj_in, sobel_y, kernel_size);
	ImageGradient gradient = get_image_gradient(G_x, G_y, width, height, kernel_size / 2);

	return gradient;
}	

void apply_gaussian_filter(GimpDrawable* drawable, gdouble sigma) {
	guint width = drawable->width;
	guint height = drawable->height;
	gint bpp = drawable->bpp;
	gint32 id = drawable->drawable_id;
		
	GimpPixelRgn rgn_in, rgn_out;
	
	gimp_pixel_rgn_init(&rgn_in, drawable, 0, 0, width, height, FALSE, FALSE);
	gimp_pixel_rgn_init(&rgn_out, drawable, 0, 0, width, height, TRUE, TRUE);
	
	gboolean has_alpha = gimp_drawable_has_alpha(id);

	guint8 kernel_size = calculate_gaussian_kernel_size(sigma);
	gdouble* kernel = g_new(gdouble, kernel_size * kernel_size);	
	calculate_gaussian_kernel(sigma, kernel_size, kernel);

	AdjustedPixelRgn* adj_in = adjust_pixel_rgn(&rgn_in, kernel_size, has_alpha);
	AdjustedPixelRgn* adj_out = get_convolution_result(adj_in, kernel, kernel_size);
	inverse_adjust_pixel_rgn(adj_out, &rgn_out, kernel_size);
	
	free_rgn(adj_out);
	free_rgn(adj_in);
	g_free(kernel);
	
	gimp_drawable_flush(drawable);
	gimp_drawable_merge_shadow(id, TRUE);
	gimp_drawable_update(id, 0, 0, width, height);
}

typedef struct {
	gdouble sigma;
	guint8 high_threshold;
	guint8 low_threshold;
} InputParams;

void run_core(GimpDrawable* drawable, gint32 image_id, const InputParams* params) {
	guint width = drawable->width;
	guint height = drawable->height;
	gint bpp = drawable->bpp;
	gint32 id = drawable->drawable_id;
			
	GimpPixelRgn rgn_out;
	gimp_pixel_rgn_init(&rgn_out, drawable, 0, 0, width, height, TRUE, TRUE);

	gimp_image_undo_enable(image_id);
	gimp_image_undo_group_start(image_id);
	
	convert_grayscale(drawable);
	apply_gaussian_filter(drawable, params->sigma);
	ImageGradient gradient = apply_sobel(drawable);
	AdjustedPixelRgn* rgn =	apply_non_max_suppression(&gradient);
	apply_hysteresis_tracking(rgn, params->high_threshold, params->low_threshold);

	inverse_adjust_pixel_rgn(rgn, &rgn_out, 0);

	free_rgn(gradient.magnitude);
	free_rgn(gradient.angle);
	free_rgn(rgn);

	gimp_drawable_flush(drawable);
	gimp_drawable_merge_shadow(id, TRUE);
	gimp_drawable_update(id, 0, 0, width, height);
	
	gimp_image_undo_group_end(image_id);	
}

static void
query (void)
{
	static GimpParamDef args[] = {
	    {
	      GIMP_PDB_INT32,
	      "run-mode",
	      "Run mode"
	    },
	    {
	      GIMP_PDB_IMAGE,
	      "image",
	      "Input image"
	    },
	    {
	      GIMP_PDB_DRAWABLE,
	      "drawable",
	      "Input drawable"
	    },
	    {
	    	GIMP_PDB_FLOAT,
	    	"sigma",
	    	"Gaussian Sigma value"
	    },
	    {
	    	GIMP_PDB_INT8,
	    	"high-threshold",
	    	"High threshold for edge detecting"
	    },
	    {
	    	GIMP_PDB_INT8,
	    	"low-threshold",
	    	"Low threshold for edge detecting"
	    }
	  };

	gimp_install_procedure (
	    "canny-edge-detector",
	    "Canny edge detector",
	    "Uses Canny method to detect edges",
	    "Saitgaliev Raif",
	    "Copyright Saitgaliev Raif",
	    "2024",
	    "_Canny",
	    "RGB*, GRAY*",
	    GIMP_PLUGIN,
	    G_N_ELEMENTS (args), 0,
	    args, NULL);

	gimp_plugin_menu_register ("canny-edge-detector",
	                               "<Image>/Filters/Misc"); 
}

static gboolean
run_plugin_dialog(GimpDrawable* drawable, InputParams* params) {
  GtkWidget *dialog;
  GtkWidget *main_vbox;
  GtkWidget *slider;
  GtkWidget *slider_label;
  GtkWidget *spin_threshold1;
  GtkWidget *spin_threshold2;
  GtkWidget *frame_label;
  gboolean   run;
  GtkWidget *entry;


  gimp_ui_init("canny_edge_detector", FALSE);

  dialog = gimp_dialog_new("Canny edge detector", "canny_edge_detector",
                           NULL, 0,
                           gimp_standard_help_func, "canny-edge_detector",
                           GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                           GTK_STOCK_OK, GTK_RESPONSE_OK, NULL);

  main_vbox = gtk_vbox_new(FALSE, 6);
  gtk_container_add(GTK_CONTAINER(GTK_DIALOG(dialog)->vbox), main_vbox);
  gtk_widget_show(main_vbox);
  
  entry = gtk_entry_new();
  gtk_widget_show(entry);
  gtk_box_pack_start(GTK_BOX(main_vbox), entry, FALSE, FALSE, 6);
  gtk_entry_set_text(GTK_ENTRY(entry), "1,0"); 

  frame_label = gtk_label_new("<b>Blurring power coefficient</b>");
  gtk_widget_show(frame_label);
  gtk_label_set_use_markup(GTK_LABEL(frame_label), TRUE);
  gtk_box_pack_start(GTK_BOX(main_vbox), frame_label, FALSE, FALSE, 6);

  // Threshold 1 Spin Button
  spin_threshold1 = gtk_spin_button_new_with_range(0, 255, 1);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin_threshold1), 100);
  gtk_widget_show(spin_threshold1);
  gtk_box_pack_start(GTK_BOX(main_vbox), spin_threshold1, FALSE, FALSE, 6);

  //Threshold 1 Label
  frame_label = gtk_label_new("<b>Threshold for sure-to-be edges (0-255)</b>"); //Corrected: Declare here
  gtk_widget_show(frame_label);
  gtk_label_set_use_markup(GTK_LABEL(frame_label), TRUE);
  gtk_box_pack_start(GTK_BOX(main_vbox), frame_label, FALSE, FALSE, 6);


  // Threshold 2 Spin Button
  spin_threshold2 = gtk_spin_button_new_with_range(0, 255, 1);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin_threshold2), 200);
  gtk_widget_show(spin_threshold2);
  gtk_box_pack_start(GTK_BOX(main_vbox), spin_threshold2, FALSE, FALSE, 6);

  //Threshold 2 Label
  frame_label = gtk_label_new("<b>Threshold for sure-to-be non-edges (0-255)</b>");
  gtk_widget_show(frame_label);
  gtk_label_set_use_markup(GTK_LABEL(frame_label), TRUE);
  gtk_box_pack_start(GTK_BOX(main_vbox), frame_label, FALSE, FALSE, 6);


  gtk_widget_show(dialog);

  run = (gimp_dialog_run(GIMP_DIALOG(dialog)) == GTK_RESPONSE_OK);

  if (run) {
    const char* sigma_str = gtk_entry_get_text(GTK_ENTRY(entry));
    params->sigma = atof(sigma_str);
    params->high_threshold = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(spin_threshold1));
    params->low_threshold = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(spin_threshold2));
  }

  gtk_widget_destroy(dialog);

  return run;
}

GimpPDBStatusType handle_input_parameters(GimpDrawable* drawable,
	GimpRunMode run_mode, 
	InputParams* params, 
	gint nparams,
	const GimpParam  *param) 
{
	GimpPDBStatusType status = GIMP_PDB_SUCCESS;
	switch(run_mode) {
		case GIMP_RUN_INTERACTIVE: {
			gimp_get_data("canny-edge-detector", params);
			if(!run_plugin_dialog(drawable, params))
				status = GIMP_PDB_CALLING_ERROR;
			break;	
		}

		case GIMP_RUN_NONINTERACTIVE: {
			if(nparams != 6) 
				status = GIMP_PDB_CALLING_ERROR;
			if(status == GIMP_PDB_SUCCESS) {
				params->sigma = param[3].data.d_float;
				params->high_threshold = param[4].data.d_float;
				params->low_threshold = param[5].data.d_float;
			}
			break;
		}
		
		case GIMP_RUN_WITH_LAST_VALS: {
			gimp_get_data("canny-edge-detector", params);
			break;
		}

		default:
			status = GIMP_PDB_CALLING_ERROR;
			break;
	}
	return status;
}

static void
run (const gchar      *name,
   gint              nparams,
   const GimpParam  *param,
   gint             *nreturn_vals,
   GimpParam       **return_vals)
{
	static GimpParam  values[1];
	GimpPDBStatusType status = GIMP_PDB_SUCCESS;
	GimpRunMode       run_mode;

	/* Setting mandatory output values */
	*nreturn_vals = 1;
	*return_vals  = values;

	run_mode = param[0].data.d_int32;

	values[0].type = GIMP_PDB_STATUS;
	values[0].data.d_status = status;
	
	GimpDrawable* drawable = gimp_drawable_get(param[2].data.d_drawable);
	gint32 image_id = param[1].data.d_image;

	InputParams params;

	status = handle_input_parameters(drawable, run_mode, &params, nparams, param);
	if(status == GIMP_PDB_SUCCESS) {
		run_core(drawable, image_id, &params);
		
		gimp_displays_flush();
		gimp_drawable_detach(drawable);


		if (run_mode == GIMP_RUN_INTERACTIVE)
		    gimp_set_data ("canny-edge-detector", &params, sizeof (InputParams));
	}	
	  return;
}

GimpPlugInInfo PLUG_IN_INFO = {
	NULL,
	NULL,
	query,
	run
};

MAIN()
