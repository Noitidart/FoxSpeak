'use strict';

// Imports
importScripts('resource://gre/modules/osfile.jsm');
importScripts('resource://gre/modules/workers/require.js');

// Globals
var core = { // have to set up the main keys that you want when aCore is merged from mainthread in init
	addon: {
		path: {
			content: 'chrome://foxspeak/content/',
		}
	},
	os: {
		name: OS.Constants.Sys.Name.toLowerCase()
	},
	firefox: {}
};

var OSStuff = {}; // global vars populated by init, based on OS

// Imports that use stuff defined in chrome
// I don't import ostypes_*.jsm yet as I want to init core first, as they use core stuff like core.os.isWinXP etc
// imported scripts have access to global vars on MainWorker.js
importScripts(core.addon.path.content + 'modules/cutils.jsm');
importScripts(core.addon.path.content + 'modules/ctypes_math.jsm');

// Setup PromiseWorker
var PromiseWorker = require('resource://gre/modules/workers/PromiseWorker.js');

// Instantiate AbstractWorker (see below).
var worker = new PromiseWorker.AbstractWorker()

worker.dispatch = function(method, args = []) {
  // Dispatch a call to method `method` with args `args`
  return self[method](...args);
};
worker.postMessage = function(...args) {
  // Post a message to the main thread
  self.postMessage(...args);
};
worker.close = function() {
  // Close the worker
  self.close();
};
worker.log = function(...args) {
  // Log (or discard) messages (optional)
  dump('Worker: ' + args.join(' ') + '\n');
};

// Connect it to message port.
self.addEventListener('message', msg => worker.handleMessage(msg));

// Define a custom error prototype.
function MainWorkerError(msgObj) {
  this.message = msgObj.message;
  this.name = msgObj.name;
}
MainWorkerError.prototype.toMsg = function() {
  return {
    exn: 'MainWorkerError',
    message: this.message,
	name: this.name
  };
};

////// end of imports and definitions

function init(objCore) { // function name init required for SIPWorker
	//console.log('in worker init');
	
	// merge objCore into core
	// core and objCore is object with main keys, the sub props
	
	core = objCore;
	
	// I import ostypes_*.jsm in init as they may use things like core.os.isWinXp etc
	switch (core.os.toolkit.indexOf('gtk') == 0 ? 'gtk' : core.os.name) {
		case 'winnt':
		case 'winmo':
		case 'wince':
			importScripts(core.addon.path.content + 'modules/ostypes_win.jsm');
			break
		case 'gtk':
			importScripts(core.addon.path.content + 'modules/ostypes_gtk.jsm');
			break;
		case 'darwin':
			importScripts(core.addon.path.content + 'modules/ostypes_mac.jsm');
			break;
		default:
			throw new MainWorkerError({
				name: 'addon-error',
				message: 'Operating system, "' + OS.Constants.Sys.Name + '" is not supported'
			});
	}
	
	// OS Specific Init
	switch (core.os.name) {
		// case 'winnt':
		// case 'winmo':
		// case 'wince':
		// 		
		// 		OSStuff.hiiii = true;
		// 		
		// 	break;
		default:
			// do nothing special
	}
	
	
	// console.log('add:', ostypes.API('add'));
	// var rez = ostypes.API('add')(5, 4);
	// console.log('rez:', rez);
    // 
	// console.log('acmod_free:', ostypes.API('acmod_free'));
	// ostypes.API('acmod_free')(null);
	
	console.log('FoxSpeak MainWorker init success');
	return true; // required for SIPWorker
}

// Start - Addon Functionality
function check_wav_header(header, expected_sr) {
    // returns 0 if invalid
    // returns 1 if valid

    if (header[34] != 0x10) {
        throw new Error('Input audio file has [%d] bits per sample instead of 16\n', header[34]);
        return 0;
    }
    if (header[20] != 0x1) {
        throw new Error('Input audio file has compression [%d] and not required PCM\n', header[20]);
        return 0;
    }
    if (header[22] != 0x1) {
        throw new Error('Input audio file has [%d] channels, expected single channel mono\n', header[22]);
        return 0;
    }
    var sr = ((header[24] & 0xFF) | ((header[25] & 0xFF) << 8) | ((header[26] & 0xFF) << 16) | ((header[27] & 0xFF) << 24));
    if (sr != expected_sr) {
        throw new Error('Input audio file has sample rate ' + sr + ', but decoder expects ' + expected_sr);
        return 0;
    }
    return 1;
}

function recognize_from_file(filename, srate) {
    // Continuous recognition from a file
    // returns undefined

    var adbuf = ostypes.TYPE.int16_t.array(2048); // int16 adbuf[2048];
    var fname = ostypes.TYPE.char(); // const char *fname;
    var hyp = ostypes.TYPE.char(); // const char *hyp;
    var k; // int32
    var utt_started = ostypes.TYPE.uint8_t(); // uint8 utt_started, in_speech;
    var in_speech = ostypes.TYPE.uint8_t(); // uint8 utt_started, in_speech;

    fname = filename;
    if ((ostypes.API('fopen')(fname, 'rb')) == null) {
        throw new Error('Failed to open file ' + fname + ' for reading');
    }

    if (fname.length > 4 && fname.substr(fname.length - 4) == '.wav') {
        var waveheader = ostypes.char.array(44); // char waveheader[44];
        ostypes.API('fread')(waveheader, 1, 44, rawfd);
        if (!check_wav_header(waveheader, srate))
            throw new Error('Failed to process file ' + fname + ' due to format mismatch.');
    }

    if (strlen(fname) > 4 && fname.substr(fname.length - 4) == '.mp3') {
        throw new Error('Can not decode mp3 files, convert input file to WAV 16kHz 16-bit mono before decoding.\n');
    }

    ostypes.API('ps_start_utt')(ps);
    utt_started = false;

    while ((k = ostypes.API('fread')(adbuf, adbuf.elementType.size, 2048, rawfd)) > 0) {
        ostypes.API('ps_process_raw')(ps, adbuf, k, false, false);
        in_speech = ostypes.API('ps_get_in_speech')(ps);
        if (in_speech && !utt_started) {
            utt_started = true;
        }
        if (!in_speech && utt_started) {
            ostypes.API('ps_end_utt')(ps);
            hyp = ostypes.API('ps_get_hyp')(ps, null, null);
            if (hyp != null) {
                console.log('hyp:', hyp);
			}
			
            ostypes.API('ps_start_utt')(ps);
            utt_started = false;
        }
    }
    ostypes.API('ps_end_utt')(ps);
    if (utt_started) {
        hyp = ostypes.API('ps_get_hyp')(ps, null, null);
        if (hyp != null) {
            console.log('hyp:', hyp);
        }
    }

    ostypes.API('fclose')(rawfd);

}

function sleep_msec(int32_ms) {
    // Sleep for specified msec
    // returns undefined

    if (core.os.name == 'winnt') {
        ostypes.API('Sleep')(ms);
    } else {
        // Unix
        var tmo = ostypes.TYPE.timeval();

        tmo.tv_sec = 0;
        tmo.tv_usec = ms * 1000;

        ostypes.API('select')(0, null, null, null, tmo.address());
    }
}
// End - Addon Functionality