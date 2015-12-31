var EXPORTED_SYMBOLS = ['ostypes'];

// no need to define core or import cutils as all the globals of the worker who importScripts'ed it are availble here

if (ctypes.voidptr_t.size == 4 /* 32-bit */) {
	var is64bit = false;
} else if (ctypes.voidptr_t.size == 8 /* 64-bit */) {
	var is64bit = true;
} else {
	throw new Error('huh??? not 32 or 64 bit?!?!');
}

var unixTypes = function() {

	// ABIs
	this.CALLBACK_ABI = ctypes.default_abi;
	this.ABI = ctypes.default_abi;
  
	// C TYPES
	this.char = ctypes.char;
	this.int = ctypes.int;
	this.int16 = ctypes.int16_t;
	this.int32 = ctypes.int32_t;
	this.long = ctypes.long;
	this.size_t = ctypes.size_t;
	this.void = ctypes.void_t;
	
	///////////// libc types
	// SIMPLE TYPES
	this.fd_set = ctypes.uint8_t; // This is supposed to be fd_set*, but on Linux at least fd_set is just an array of bitfields that we handle manually. link4765403
	
	// SUPER ADVANCED TYPES // defined by "advanced types"

	// SUPER DUPER ADVANCED TYPES // defined by "super advanced types"

	// STRUCTURES
	// consts for structures
	var struct_const = {
		
	};
	
	// SIMPLE STRUCTS // based on any of the types above
	this.FILE = ctypes.StructType('FILE');
	this.timeval = ctypes.StructType('timeval', [
		{ 'tv_sec': this.long },
		{ 'tv_usec': this.long }
	]);
	
	// ADVANCED STRUCTS // based on "simple structs" to be defined first

	// FURTHER ADVANCED STRUCTS
	
	// FURTHER ADV STRUCTS

	// FUNCTION TYPES
	
	
	// STRUCTS USING FUNC TYPES
	
	
	
	
	/////// pocket-sphinx types
	this.acmod_t = ctypes.void_t;
	this.ps_decoder_t = ctypes.StructType('ps_decoder_t');
	
	// this.ps_decoder_t = ctypes.StructType('ps_decoder_t', [ // http://cmusphinx.sourceforge.net/doc/pocketsphinx/structps__decoder__s.html
	// 	{ config: this.cmd_ln_t.ptr },
	// 	{ refcount: this.int },
	// 	{ acmod: this.acmod_t.ptr },
	// 	{ dict: this.dict_t.ptr },
	// 	{ d2p: this.dict2pid_t.ptr },
	// 	{ lmath: this.logmath_t.ptr },
	// 	{ searches: this.hash_table_t.ptr },
	// 	{ search: this.ps_search_t.ptr },
	// 	{ phone_loop: this.ps_search_t.ptr },
	// 	{ pl_window: this.int },
	// 	{ uttno: this.uint32 },
	// 	{ perf: this.ptmr_t },
	// 	{ n_frame: this.uint32 },
	// 	{ mfclogdir: this.char.ptr },
	// 	{ rawlogdir: this.char.ptr },
	// 	{ senlogdir: this.char.ptr }
	// ]);

};

var unixInit = function() {
	var self = this;

	this.IS64BIT = is64bit;

	this.TYPE = new unixTypes();

	// CONSTANTS
	this.CONST = {
		
	};

	var _lib = {}; // cache for lib
	var lib = function(path) {
		//ensures path is in lib, if its in lib then its open, if its not then it adds it to lib and opens it. returns lib
		//path is path to open library
		//returns lib so can use straight away

		if (!(path in _lib)) {
			//need to open the library
			//default it opens the path, but some things are special like libc in mac is different then linux or like x11 needs to be located based on linux version
			switch (path) {
				case 'add':
					
						var nativeFileExtension;
						if (core.os.name == 'darwin') {
							nativeFileExtension = 'dylib';
						} else if (core.os.toolkit.indexOf('gtk') == 0) {
							nativeFileExtension = 'so';
						} else {
							nativeFileExtension = 'dll';
						}
						console.log('will now try to open add', OS.Path.join(core.addon.path.file, 'modules', 'add.' + nativeFileExtension));
						_lib[path] = ctypes.open(OS.Path.join(core.addon.path.file, 'modules', 'add.' + nativeFileExtension));

					
					break;
				case 'pocket-sphinx':
					
						var nativeFileExtension;
						if (core.os.name == 'darwin') {
							nativeFileExtension = 'dylib';
						} else if (core.os.toolkit.indexOf('gtk') == 0) {
							nativeFileExtension = 'so';
						} else {
							nativeFileExtension = 'dll';
						}
						
						_lib[path] = ctypes.open(OS.Path.join(core.addon.path.file, 'modules', 'libpocketsphinx.' + nativeFileExtension));
						
					
					break;
				case 'libc':

						switch (core.os.name) {
							case 'darwin':
								_lib[path] = ctypes.open('libc.dylib');
								break;
							case 'freebsd':
								_lib[path] = ctypes.open('libc.so.7');
								break;
							case 'openbsd':
								_lib[path] = ctypes.open('libc.so.61.0');
								break;
							case 'android':
							case 'sunos':
							case 'netbsd': // physically unverified
							case 'dragonfly': // physcially unverified
								_lib[path] = ctypes.open('libc.so');
								break;
							case 'linux':
								_lib[path] = ctypes.open('libc.so.6');
								break;
							case 'gnu/kfreebsd': // physically unverified
								_lib[path] = ctypes.open('libc.so.0.1');
								break;
							default:
								throw new MainWorkerError({
									name: 'addon-error',
									message: 'Path to libc on operating system of , "' + OS.Constants.Sys.Name + '" is not supported'
								});
						}

					break;
				default:
					try {
						_lib[path] = ctypes.open(path);
					} catch (ex) {
						throw new Error({
							name: 'addon-error',
							message: 'Could not open ctypes library path of "' + path + '"',
							ex_msg: ex.message
						});
					}
			}
		}
		return _lib[path];
	};

	// start - function declares
	var _api = {};
	this.API = function(declaration) { // it means ensureDeclared and return declare. if its not declared it declares it. else it returns the previously declared.
		if (!(declaration in _api)) {
			_api[declaration] = preDec[declaration](); //if declaration is not in preDec then dev messed up
		}
		return _api[declaration];
	};

	// start - predefine your declares here
	var preDec = { //stands for pre-declare (so its just lazy stuff) //this must be pre-populated by dev // do it alphabateized by key so its ez to look through
		// start - libc
		fclose: function() {
			/* http://linux.die.net/man/3/fclose
			 * int fclose(
			 *   FILE *fp
			 * );
			 */
			return lib('libc').declare('fclose', self.TYPE.ABI,
				self.TYPE.int,		// return
				self.TYPE.FILE.ptr	// *fp
			)
		},
		fopen: function() {
			/* http://man7.org/linux/man-pages/man3/fopen.3.html
			 * FILE *fopen(
			 *   const char *path,
			 *   const char *mode
			 * );
			 */
			return lib('libc').declare('fopen', self.TYPE.ABI,
				self.TYPE.FILE.ptr,		// return
				self.TYPE.char.ptr,		// *path
				self.TYPE.char.ptr		// *mode
			);
		},
		fread: function() {
			/* http://linux.die.net/man/3/fread
			 * size_t fread (
			 *   void *ptr,
			 *   size_t size,
			 *   size_t nmemb,
			 *   FILE *stream
			 * );
			 */
			return lib('libc').declare('fread', self.TYPE.ABI, 
				self.TYPE.size_t,		// return
				self.TYPE.void.ptr,		// *ptr
				self.TYPE.size_t, 		// size
				self.TYPE.size_t, 		// count
				self.TYPE.FILE.ptr		// *stream
			);
		},
		select: function() {
			/* http://linux.die.net/man/2/select
			 * int select (
			 *   int nfds,
			 *   fd_set *readfds,
			 *   fd_set *writefds,
			 *   fd_set *exceptfds,
			 *   struct timeval *timeout
			 * );
			 */
			return lib('libc').declare('select', self.TYPE.ABI,
				self.TYPE.int,			// return
				self.TYPE.int,			// nfds
				self.TYPE.fd_set.ptr,	// *readfds  // This is supposed to be fd_set*, but on Linux at least fd_set is just an array of bitfields that we handle manually. link4765403
				self.TYPE.fd_set.ptr,	// *writefds // This is supposed to be fd_set*, but on Linux at least fd_set is just an array of bitfields that we handle manually. link4765403
				self.TYPE.fd_set.ptr,	// *exceptfds // This is supposed to be fd_set*, but on Linux at least fd_set is just an array of bitfields that we handle manually. link4765403
				self.TYPE.timeval.ptr	// *timeout
			);
		},
		// end - libc
		// start - pocket-sphinx
		acmod_free: function() {
			/* 
			 * void acmod_free (
			 *   acmod_t *acmod
			 * );
			 */
			return lib('pocket-sphinx').declare('acmod_free', self.TYPE.ABI,
				self.TYPE.void,			// return
				self.TYPE.acmod_t.ptr	// *acmod_t
			);
		},
		ps_end_utt: function() {
			/* 
			 * int ps_end_utt(ps_decoder_t *ps)
			 */
			return lib('pocket-sphinx').declare('ps_end_utt', self.TYPE.ABI,
				self.TYPE.int,				// return
				self.TYPE.ps_decoder_t.ptr	// *ps
			);
		},
		ps_get_hyp: function() {
			/*
			 * char const * ps_get_hyp(ps_decoder_t *ps, int32 *out_best_score, char const **out_uttid)
			 */
			return lib('pocket-sphinx').declare('ps_get_hyp', self.TYPE.ABI,
				self.TYPE.char.ptr,			// return
				self.TYPE.ps_decoder_t.ptr,	// *ps
				self.TYPE.int32.ptr,		// *out_best_score
				self.TYPE.char.ptr.ptr		// **out_uttid
			);
		},
		ps_process_raw: function() {
			/* 
			 * int ps_process_raw(ps_decoder_t *ps,
             *       int16 const *data,
             *       size_t n_samples,
             *       int no_search,
             *       int full_utt);
			 */
			return lib('pocket-sphinx').declare('ps_process_raw', self.TYPE.ABI,
				self.TYPE.int,			// return
				self.TYPE.int16.ptr,	// *data
				self.TYPE.size_t,		// n_samples
				self.TYPE.int,			// no_search
				self.TYPE.int			// full_utt
			);
		},
		ps_start_utt: function() {
			/* 
			 * int ps_start_utt(ps_decoder_t *ps, char const *uttid)
			 */
			return lib('pocket-sphinx').declare('ps_start_utt', self.TYPE.ABI,
				self.TYPE.int,				// return
				self.TYPE.ps_decoder_t.ptr,	// *ps
				self.TYPE.char.ptr			// *uttid
			);
		},
		// end - pocket-sphinx
		add: function() {
			return lib('add').declare('add', self.TYPE.ABI,
				ctypes.int32_t, // return type
				ctypes.int32_t, // arg1 type
				ctypes.int32_t // arg2 type
			);
		}
	};
	// end - predefine your declares here
	// end - function declares

	this.HELPER = {

	};
};

var ostypes = new unixInit();