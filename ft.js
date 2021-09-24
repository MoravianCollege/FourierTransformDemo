/* global rfft2d, irfft2d */
/* exported $, setup_canvas, set_canvas_data, get_image_bytes, download_image */
/* exported sine2d, solid_gray, to_grayscale, from_grayscale, average_waves_and_fft, fft_to_image */
/* exported compute_fft, to_f32, get_fft_info, compute_sorted_fft_indices */
"use strict";

const WIDTH = 128, HEIGHT = 128, SCALE = 2;


/** Convience function for getting an element by id. */
function $(id) { return document.getElementById(id); }


/** Find a canvas by ID, set its size, and get the 2D context. */
function setup_canvas(id) {
	let canvas = document.getElementById(id);
	canvas.width = WIDTH;
	canvas.height = HEIGHT;
	let context = canvas.getContext('2d');
	context.imageSmoothingEnabled = false;
	return context;
}


/** Set a canvas (from its context) to display the given image. */
function set_canvas_data(context, im) {
	context.clearRect(0, 0, context.canvas.width, context.canvas.height);
	//context.canvas.width += 0;
	if (!('im_data' in context)) {
		context.im_data = context.createImageData(context.canvas.width, context.canvas.height);
	}
	context.im_data.data.set(im);
	context.putImageData(context.im_data, 0, 0);
}


/**
 * Get the RGBA bytes from an <img> or Image object. The image is forced to be
 * width-by-height.
 */
 function get_image_bytes(img) {
	img.width = WIDTH;
	img.height = HEIGHT;
	let canvas = document.createElement('canvas');
	canvas.width = WIDTH;
	canvas.height = HEIGHT;
	let context = canvas.getContext('2d');
	context.drawImage(img, 0, 0, WIDTH, HEIGHT);
	return context.getImageData(0, 0, WIDTH, HEIGHT).data;
}


/**
 * Download the image being shown in the given canvas as a PNG file.
 * The canvas can either be a <canvas> or the context.
 */
 function download_image(canvas) {
	let a = document.createElement('a'); // a dummy <a> element
	if ('canvas' in canvas) { canvas = canvas.canvas; }
	a.href = canvas.toDataURL(); // creates a PNG image
	a.setAttribute('download', 'image.png');
	a.click();
}


/**
 * Create a 2D sine wave in an image buffer with the given properties:
 *   - amplitude, in gray levels (from 0 to 255 is reasonable)
 *   - frequency, in cycles per pixel/sample (from 0 to 0.5 is reasonable)
 *   - angle, in radians (from -pi to +pi is reasonable)
 *   - phase, in radians (from -pi to +pi is reasonable, default 0)
 *   - offset, in gray levels (from 0 to 255 is reasonable, default is 127.5)
 * This uses the global width and height properties for the image size.
 * This returns a 1D array of uint8 values.
 */
function sine2d(amp, freq, angle, phase = 0, offset = 127.5) {
	let im = new Uint8ClampedArray(WIDTH * HEIGHT * 4);
	let rad_x = 2*Math.PI*freq*Math.cos(angle);
	let rad_y = 2*Math.PI*freq*Math.sin(angle);
	for (let y = 0; y < HEIGHT; y++) {
		for (let x = 0; x < WIDTH; x++) {
			let pos = (y * WIDTH + x) * 4;
			let value = Math.round(amp * Math.sin(rad_x*x+rad_y*y+phase) + offset);
			im[pos  ] = value;
			im[pos+1] = value;
			im[pos+2] = value;
			im[pos+3] = 255;
		}
	}
	return im;
}


/**
 * Generates a solid-color gray image.
 */
function solid_gray(value) {
	let im = new Uint8ClampedArray(WIDTH * HEIGHT * 4);
	for (let y = 0; y < HEIGHT; y++) {
		for (let x = 0; x < WIDTH; x++) {
			let pos = (y * WIDTH + x) * 4;
			im[pos  ] = value;
			im[pos+1] = value;
			im[pos+2] = value;
			im[pos+3] = 255;
		}
	}
	return im;
}


/**
 * Average the array of waves into a single image and perform an FFT.
 * This returns the averaged image and the FFT data.
 */
function average_waves_and_fft(waves) {
	let im = new Uint8ClampedArray(WIDTH * HEIGHT * 4);
	let im_f32 = new Float32Array(WIDTH * HEIGHT); // need raw gray data for FFT
	if (waves.length !== 0) {
		for (let y = 0; y < HEIGHT; y++) {
			for (let x = 0; x < WIDTH; x++) {
				let pos_raw = y * WIDTH + x;
				let pos = pos_raw * 4;
				let value = 0;
				for (let i = 0; i < waves.length; i++) { value += waves[i][pos]; }
				value /= waves.length;
				im_f32[pos_raw] = value;
				im[pos] = im[pos+1] = im[pos+2] = Math.round(value);
				im[pos+3] = 255;
			}
		}
	}
	return [im, rfft2d(im_f32, HEIGHT, WIDTH)];
}


/** Convert an RGBA image to grayscale using Rec 709. */
function to_grayscale(im) {
	let gray = new Uint8ClampedArray(WIDTH * HEIGHT);
	for (let y = 0; y < HEIGHT; y++) {
		for (let x = 0; x < WIDTH; x++) {
			let pos = y*WIDTH + x, pos4 = pos*4;
			gray[pos] = (im[pos4]*0.2126 + im[pos4+1]*0.7152 + im[pos4+2]*0.0722) * im[pos4+3]/255;
		}
	}
	return gray;
}


/** Convert a grayscale image to RGBA. */
function from_grayscale(im) {
	let rgb = new Uint8ClampedArray(WIDTH * HEIGHT * 4);
	for (let y = 0; y < HEIGHT; y++) {
		for (let x = 0; x < WIDTH; x++) {
			let pos = y*WIDTH + x, pos4 = pos*4;
			rgb[pos4] = rgb[pos4+1] = rgb[pos4+2] = im[pos];
			rgb[pos4+3] = 255;
		}
	}
	return rgb;
}


/** Convert array to float-32. */
function to_f32(data) {
	let out = new Float32Array(data.length);
	for (let i = 0; i < data.length; i++) { out[i] = data[i]; }
	return out;
}


/** Compute the FFT of the image. */
function compute_fft(im) {
	if (im.constructor.name !== 'Float32Array') {
		im = to_f32(im);
	}
	return rfft2d(im, HEIGHT, WIDTH);
}


/**
 * Takes an FFT spectrum that contains complex data which alternates between
 * real and imaginary values and converts it into an image. The image is
 * colored using the HSV (hue, saturation, value) model with the value being
 * based on the log of the magnitude of the complex values, the hue based on
 * the angle of the complex number, and the saturation is always 1.
 * 
 * This also returns the minimum magnitude and the maximum magnitude (both
 * un-scaled and un-logged).
 */
function fft_to_image(fft, color = false, log = true) {
	const W2 = WIDTH/2 + 1;
	let n = HEIGHT*W2; // TODO: fft.length/2; // real and imaginary interleaved;
	let mag = new Float32Array(n); // magnitude (or log-of-magnitude, not normalized)
	let hue = new Float32Array(n); // hue (derived from angle)
	let max = 0, min = 1e16; // max/min (not including DC component)
	for (let i = 1; i < n; i++) {
		let real = fft[2*i], imag = fft[2*i+1];
		let val = Math.hypot(real, imag);
		if (log) { val = Math.log1p(val); }
		if (val > max) { max = val; } else if (val < min) { min = val; }
		mag[i] = val;
		if (color) {
			let h = Math.atan2(imag, real); //complex angle in the range [-pi, pi)
			if (h < 0) { h += 2*Math.PI; } // now in the range [0, 2*pi)
			if (!isFinite(h)) { h = 0; } // just in case
			hue[i] = (h/(2*Math.PI)) % 1; // now in the range [0, 1)
		}
	}
	let dc = Math.hypot(fft[0], fft[1]);
	dc = log ? Math.log1p(dc + 1) : dc;
	mag[0] = dc > max ? max : dc; // clip DC component
	hue[0] = 0; // correct DC component angle

	const scale = 255/(max-min), size = WIDTH*HEIGHT;
	let im = new Uint8ClampedArray(size * 4);
	for (let i = 0; i < n; i++) {
		// Get the pixel color
		let val = (mag[i]-min)*scale;
		let r1, g1, b1, r2, g2, b2;
		if (color && i !== 0) {
			[r1, g1, b1] = h1v_to_rgb(hue[i], val);
			[r2, g2, b2] = h1v_to_rgb(1-hue[i], val);
		} else {
			[r1, g1, b1] = [val, val, val];
			[r2, g2, b2] = [val, val, val];
		}

		// Get the position of the pixel
		// We have to apply an fft-shift along with dealing with the symmetries in an rfft
		// TODO: for some reason half of the first row is missing...
		let x = i % W2 + WIDTH/2, y = (Math.trunc(i/W2)+HEIGHT/2)%HEIGHT;
		let j = y*WIDTH+x, k = size + WIDTH - j;

		// Set the pixels
		im[j*4  ] = Math.round(r1);
		im[j*4+1] = Math.round(g1);
		im[j*4+2] = Math.round(b1);
		im[k*4  ] = Math.round(r2);
		im[k*4+1] = Math.round(g2);
		im[k*4+2] = Math.round(b2);
		im[j*4+3] = im[k*4+3] = 255; // alpha
	}

	return [im, log ? Math.expm1(min) : min, log ? Math.expm1(max) : max];
}


/** Convert HSV color to RGB color but assuming S = 1 */
function h1v_to_rgb(h, v) {
	let i = Math.trunc(h*6.0);
	let f = h*6.0 - i;
	let t = v*(i%2 === 0 ? f : 1.0-f);
	i %= 6;
	if	  (i === 0) return [v, t, 0];
	else if (i === 1) return [t, v, 0];
	else if (i === 2) return [0, v, t];
	else if (i === 3) return [0, t, v];
	else if (i === 4) return [t, 0, v];
	else /*if (i === 5)*/ return [v, 0, t];
}


/**
 * Retrieves information about a location in the FFT.
 */
function get_fft_info(fft, x, y) {
	const W2 = WIDTH/2+1;

	let freq = Math.hypot(x, y) / Math.max(WIDTH, HEIGHT);
	let angle = Math.atan2(y, x);
	let mirrored = x < 0;
	if (y < 0) { y += HEIGHT; } // change from [-h/2 to h/2) to [0 to h)
	if (mirrored) { x = -x; y = (HEIGHT-y)%HEIGHT; } // on the mirror side
	
	let dc = Math.hypot(fft[0], fft[1])/(WIDTH*HEIGHT);
	if (Math.abs(freq) < 0.005) { return [0, 0, dc, 0, solid_gray(dc)]; }

	let pos = 2*(y*W2+x%W2); // position in the FFT array, 2x since it is real/imag interleaved
	let real = fft[pos], imag = fft[pos+1];
	if (mirrored) { imag = -imag; }

	// TODO:
	let sum = 0; // everything but the DC component
	for (let i = 2; i < W2*HEIGHT*2; i+=2) { sum += Math.hypot(current_fft[i], current_fft[i+1]); }

	let value = Math.hypot(real, imag); // / sum * dc; // in range [0, sqrt(width*height)*255]
	let phase = Math.atan2(imag, real); // in range [-pi, pi)

	//let im = sine2d(value, freq, angle, phase, value);
	let partial = extract_fft(fft, [pos]);    // Perform inverse FFT
    let wave = from_grayscale(irfft2d(partial, HEIGHT, WIDTH));

	return [freq, angle, value, phase, wave];
}

/**
 * Extract a subset of an FFT.
 */
function extract_fft(fft, indices, count = indices.length) {
	let partial = new Float32Array((WIDTH/2+1)*HEIGHT*2);
	// Force the DC component
	partial[0] = fft[0];
	partial[1] = fft[1];
	// Fill in the partial FFT
	for (let i = 0; i < count; i++) {
		let pos = indices[i] * 2;
		partial[pos] = fft[pos];
		partial[pos+1] = fft[pos+1];
	}
	return partial;
}

/**
 * Gets a list of the highest amplitude indices in an FFT.
 */
function compute_sorted_fft_indices(fft) {
	let N = (WIDTH / 2 + 1) * HEIGHT;
	let indices = new Array(N);
	let magnitudes = new Float32Array(N);
	for (let i = 0; i < N; i++) {
		indices[i] = i;
		magnitudes[i] = Math.hypot(fft[2*i], fft[2*i+1]);
	}
	indices.sort((a, b) => {
		if (magnitudes[a] < magnitudes[b]) { return 1; }
		if (magnitudes[a] > magnitudes[b]) { return -1; }
		return 0;
	});
	return indices;
}


/**
 * Computes a 'histogram' of the frequency magnitudes in an FFT. The magnitudes
 * are summed for all things in that bucket, all the way around the circle.
 * Each histogram bucket represents 2 "circles" of frequencies, divide the
 * index by Math.max(WIDTH, HEIGHT)*2 to get the frequency of the bucket.
 * 
 * The reason for 2 "circles" per bucket is that due to aliasing issues, every
 * other circle will have many more points in it and this smooths that out.
 */
function freq_hist(fft) {
	let hist = new Float64Array(Math.round(Math.hypot(WIDTH, HEIGHT)/4) + 1);
	const W2 = WIDTH / 2 + 1;
	for (let y = -HEIGHT/2; y < HEIGHT/2; y++) {
		for (let x = -WIDTH/2; x <= 0; x++) {
			let x_ = x, y_ = y < 0 ? y + HEIGHT : y;
			if (x < 0) { x_ = -x; y_ = HEIGHT - y_; }
			let pos = 2*(y_*W2+x_%W2);
			let dist = Math.round(Math.hypot(x, y)/2);
			let mag = Math.hypot(fft[pos], fft[pos+1]);
			hist[dist] += mag;
			if (x != 0) { hist[dist] += mag; }
		}
	}
	return hist;
}
