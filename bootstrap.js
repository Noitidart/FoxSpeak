// Imports
const {interfaces: Ci, utils: Cu, classes:Cc} = Components;
Cu.import('resource://gre/modules/Services.jsm');
Cu.import('resource:///modules/CustomizableUI.jsm');

// Globals
var core = {
	addon: {
		name: 'FoxSpeak',
		id: 'FoxSpeak@jetpack',
		path: {
			name: 'foxspeak',
			//
			content: 'chrome://foxspeak/content/',
			locale: 'chrome://foxspeak/locale/',
			//
			modules: 'chrome://foxspeak/content/modules/',
			workers: 'chrome://foxspeak/content/modules/workers/',
			//
			resources: 'chrome://foxspeak/content/resources/',
			images: 'chrome://foxspeak/content/resources/images/',
			scripts: 'chrome://foxspeak/content/resources/scripts/',
			styles: 'chrome://foxspeak/content/resources/styles/',
			fonts: 'chrome://foxspeak/content/resources/styles/fonts/',
			pages: 'chrome://foxspeak/content/resources/pages/'
			// below are added by worker
			// storage: OS.Path.join(OS.Constants.Path.profileDir, 'jetpack', core.addon.id, 'simple-storage')
		},
		pref_branch: 'extensions.FoxSpeak@jetpack.',
		cache_key: Math.random() // set to version on release
	},
	os: {
		// // name: added by worker
		// // mname: added by worker
		toolkit: Services.appinfo.widgetToolkit.toLowerCase(),
		xpcomabi: Services.appinfo.XPCOMABI
	},
	firefox: {
		// pid: Services.appinfo.processID,
		// version: Services.appinfo.version,
		// channel: Services.prefs.getCharPref('app.update.channel')
	}
};

var gWkComm;
var gGenCssUri;
var gCuiCssUri;

var windowListener = {
	//DO NOT EDIT HERE
	onOpenWindow: function (aXULWindow) {
		// Wait for the window to finish loading
		var aDOMWindow = aXULWindow.QueryInterface(Ci.nsIInterfaceRequestor).getInterface(Ci.nsIDOMWindow);
		aDOMWindow.addEventListener('load', function () {
			aDOMWindow.removeEventListener('load', arguments.callee, false);
			windowListener.loadIntoWindow(aDOMWindow);
		}, false);
	},
	onCloseWindow: function (aXULWindow) {},
	onWindowTitleChange: function (aXULWindow, aNewTitle) {},
	register: function () {

		// Load into any existing windows
		let DOMWindows = Services.wm.getEnumerator(null);
		while (DOMWindows.hasMoreElements()) {
			let aDOMWindow = DOMWindows.getNext();
			if (aDOMWindow.document.readyState == 'complete') { //on startup `aDOMWindow.document.readyState` is `uninitialized`
				windowListener.loadIntoWindow(aDOMWindow);
			} else {
				aDOMWindow.addEventListener('load', function () {
					aDOMWindow.removeEventListener('load', arguments.callee, false);
					windowListener.loadIntoWindow(aDOMWindow);
				}, false);
			}
		}
		// Listen to new windows
		Services.wm.addListener(windowListener);
	},
	unregister: function () {
		// Unload from any existing windows
		let DOMWindows = Services.wm.getEnumerator(null);
		while (DOMWindows.hasMoreElements()) {
			let aDOMWindow = DOMWindows.getNext();
			windowListener.unloadFromWindow(aDOMWindow);
		}
		/*
		for (var u in unloaders) {
			unloaders[u]();
		}
		*/
		//Stop listening so future added windows dont get this attached
		Services.wm.removeListener(windowListener);
	},
	//END - DO NOT EDIT HERE
	loadIntoWindow: function (aDOMWindow) {
		if (!aDOMWindow) { return }

		if (aDOMWindow.gBrowser) {
			var domWinUtils = aDOMWindow.QueryInterface(Ci.nsIInterfaceRequestor).getInterface(Ci.nsIDOMWindowUtils);
			domWinUtils.loadSheet(gCuiCssUri, domWinUtils.AUTHOR_SHEET);
			// domWinUtils.loadSheet(gGenCssUri, domWinUtils.AUTHOR_SHEET);
		}
	},
	unloadFromWindow: function (aDOMWindow) {
		if (!aDOMWindow) { return }

		if (aDOMWindow.gBrowser) {
			var domWinUtils = aDOMWindow.QueryInterface(Ci.nsIInterfaceRequestor).getInterface(Ci.nsIDOMWindowUtils);
			domWinUtils.removeSheet(gCuiCssUri, domWinUtils.AUTHOR_SHEET);
			// domWinUtils.removeSheet(gGenCssUri, domWinUtils.AUTHOR_SHEET);
		}
	}
};

var gRecording;
var gRecordingDuration;
function cuiClick(e) {
	if (!gRecording) {
		var w = Services.appShell.hiddenDOMWindow;
		w.navigator.mediaDevices.getUserMedia({
			audio: true
		}).then(function (stream) {
			// do something with the stream
			console.log('ok got stream:', stream);
			gRecording = new w.MediaRecorder(stream);
			gRecording.mimeType = 'audio/ogg';

			gRecording.addEventListener('dataavailable', function(e) {
				console.log('data avail, e:')

				var fileReader = new w.FileReader();
				fileReader.onload = function() {
					var arrbuf = this.result;

					gWkComm.postMessage('readRecording', arrbuf, [arrbuf], function(url, aComm) {
						console.log('ok got url, play it');
						var w = Services.wm.getMostRecentWindow('navigator:browser');
						var audioEl = w.document.createElementNS('http://www.w3.org/1999/xhtml', 'audio');
						audioEl.setAttribute('autoplay', 'true');
						// audioEl.addEventListener('end', function() {
						// 	console.log('audio ended');
							// gWkComm.postMessage('releaseUrl', url)
							// audioEl.parentNode.removeChild(audioEl);
						// });
						// w.document.documentElement.appendChild(audioEl);
						audioEl.src = url;
						w.setTimeout(function() {
							gWkComm.postMessage('releaseUrl', url);
						}, gRecordingDuration + 1000);
					});

					gRecording = null;
				};
				fileReader.readAsArrayBuffer(e.data);
			}, false);

			gRecordingDuration = Date.now();
			gRecording.start();

			console.log('recording started');
		}, function(aReason) {
			console.error('failed, aReason:', aReason);
		});
	} else {
		console.log('recording stopped');
		gRecordingDuration = Date.now() - gRecordingDuration;
		gRecording.stop();
	}
}

// bootstrap
function install() {}

function uninstall(aData, aReason) {
	if (aReason == ADDON_UNINSTALL) {
		// delete prefs
	}
}

function startup(aData, aReason) {

	// set preferences defaults

	gWkComm = new workerComm(core.addon.path.workers + 'MainWorker.js', ()=>{return core}, function(aArg, aComm) {
		core = aArg;

		gCuiCssUri = Services.io.newURI(core.addon.path.styles + 'cui.css', null, null);
		// gGenCssUri = Services.io.newURI(core.addon.path.styles + 'general.css', null, null);

		// window listener
		windowListener.register();

		// insert cui
		CustomizableUI.createWidget({
			id: 'cui_foxspeak',
			defaultArea: CustomizableUI.AREA_NAVBAR,
			label: formatStringFromNameCore('cui_foxspeak_lbl', 'global'),
			tooltiptext: formatStringFromNameCore('cui_foxspeak_tip', 'global'),
			onCommand: cuiClick
		});

	});

	gWkComm.postMessage('dummyForInstantInstantiate');

}

function shutdown(aData, aReason) {

	if (aReason == APP_SHUTDOWN) { return }

	CustomizableUI.destroyWidget('cui_foxspeak');

	windowListener.unregister();

	workerComm_unregAll();
}

// start - functions called by worker

// end - functions called by worker

//start - common helper functions
function formatStringFromNameCore(aLocalizableStr, aLoalizedKeyInCoreAddonL10n, aReplacements) {
	// 051916 update - made it core.addon.l10n based
    // formatStringFromNameCore is formating only version of the worker version of formatStringFromName, it is based on core.addon.l10n cache

	// try {
	// 	var cLocalizedStr = core.addon.l10n[aLoalizedKeyInCoreAddonL10n];
	// } catch (ex) {
	// 	console.error('formatStringFromNameCore error:', ex, 'args:', aLocalizableStr, aLoalizedKeyInCoreAddonL10n, aReplacements);
	// }
	var cLocalizedStr = core.addon.l10n[aLoalizedKeyInCoreAddonL10n][aLocalizableStr];
	// console.log('cLocalizedStr:', cLocalizedStr, 'args:', aLocalizableStr, aLoalizedKeyInCoreAddonL10n, aReplacements);
    if (aReplacements) {
        for (var i=0; i<aReplacements.length; i++) {
            cLocalizedStr = cLocalizedStr.replace('%S', aReplacements[i]);
        }
    }

    return cLocalizedStr;
}

//rev1 - https://gist.github.com/Noitidart/c4ab4ca10ff5861c720b
function validateOptionsObj(aOptions, aOptionsDefaults) {
	// ensures no invalid keys are found in aOptions, any key found in aOptions not having a key in aOptionsDefaults causes throw new Error as invalid option
	for (var aOptKey in aOptions) {
		if (!(aOptKey in aOptionsDefaults)) {
			console.error('aOptKey of ' + aOptKey + ' is an invalid key, as it has no default value, aOptionsDefaults:', aOptionsDefaults, 'aOptions:', aOptions);
			throw new Error('aOptKey of ' + aOptKey + ' is an invalid key, as it has no default value');
		}
	}

	// if a key is not found in aOptions, but is found in aOptionsDefaults, it sets the key in aOptions to the default value
	for (var aOptKey in aOptionsDefaults) {
		if (!(aOptKey in aOptions)) {
			aOptions[aOptKey] = aOptionsDefaults[aOptKey];
		}
	}
}

function Deferred() { // revFinal
	this.resolve = null;
	this.reject = null;
	this.promise = new Promise(function(resolve, reject) {
		this.resolve = resolve;
		this.reject = reject;
	}.bind(this));
	Object.freeze(this);
}
function genericReject(aPromiseName, aPromiseToReject, aReason) {
	var rejObj = {
		name: aPromiseName,
		aReason: aReason
	};
	console.error('Rejected - ' + aPromiseName + ' - ', rejObj);
	if (aPromiseToReject) {
		aPromiseToReject.reject(rejObj);
	}
}
function genericCatch(aPromiseName, aPromiseToReject, aCaught) {
	var rejObj = {
		name: aPromiseName,
		aCaught: aCaught
	};
	console.error('Caught - ' + aPromiseName + ' - ', rejObj);
	if (aPromiseToReject) {
		aPromiseToReject.reject(rejObj);
	}
}
// start - CommAPI
// common to all of these apis
	// whenever you use the message method, the method MUST not be a number, as if it is, then it is assumed it is a callback
	// if you want to do a transfer of data from a callback, if transferring is supported by the api, then you must wrapp it in aComm.CallbackTransferReturn

var gBootstrap = this;

// start - CommAPI for bootstrap-framescript - bootstrap side - cross-file-link55565665464644
// message method - transcribeMessage - it is meant to indicate nothing can be transferred, just copied/transcribed to the other process
// first arg to transcribeMessage is a message manager, this is different from the other comm api's
var gCrossprocComms = [];
function crossprocComm_unregAll() {
	var l = gCrossprocComms.length;
	for (var i=0; i<l; i++) {
		gCrossprocComms[i].unregister();
	}
}
function crossprocComm(aChannelId) {
	// when a new framescript creates a crossprocComm on framscript side, it requests whatever it needs on init, so i dont offer a onBeforeInit or onAfterInit on bootstrap side

	var scope = gBootstrap;
	gCrossprocComms.push(this);

	this.unregister = function() {
		Services.mm.removeMessageListener(aChannelId, this.listener);

		var l = gCrossprocComms.length;
		for (var i=0; i<l; i++) {
			if (gCrossprocComms[i] == this) {
				gCrossprocComms.splice(i, 1);
				break;
			}
		}

		// kill framescripts
		Services.mm.broadcastAsyncMessage(aChannelId, {
			method: 'UNINIT_FRAMESCRIPT'
		});
	};

	this.listener = {
		receiveMessage: function(e) {
			var messageManager = e.target.messageManager;
			var browser = e.target;
			var payload = e.data;
			console.log('bootstrap crossprocComm - incoming, payload:', payload); //, 'e:', e);
			// console.log('this in receiveMessage bootstrap:', this);

			if (payload.method) {
				if (!(payload.method in scope)) { console.error('method of "' + payload.method + '" not in scope'); throw new Error('method of "' + payload.method + '" not in scope') }  // dev line remove on prod
				var rez_bs_call = scope[payload.method](payload.arg, messageManager, browser, this); // only on bootstrap side, they get extra 2 args
				if (payload.cbid) {
					if (rez_bs_call && rez_bs_call.constructor.name == 'Promise') {
						rez_bs_call.then(
							function(aVal) {
								console.log('Fullfilled - rez_bs_call - ', aVal);
								this.transcribeMessage(messageManager, payload.cbid, aVal);
							}.bind(this),
							genericReject.bind(null, 'rez_bs_call', 0)
						).catch(genericCatch.bind(null, 'rez_bs_call', 0));
					} else {
						console.log('calling transcribeMessage for callbck with args:', payload.cbid, rez_bs_call);
						this.transcribeMessage(messageManager, payload.cbid, rez_bs_call);
					}
				}
			} else if (!payload.method && payload.cbid) {
				// its a cbid
				this.callbackReceptacle[payload.cbid](payload.arg, messageManager, browser, this);
				delete this.callbackReceptacle[payload.cbid];
			} else {
				throw new Error('invalid combination');
			}
		}.bind(this)
	};
	this.nextcbid = 1; //next callback id
	this.transcribeMessage = function(aMessageManager, aMethod, aArg, aCallback) {
		// console.log('bootstrap sending message to framescript', aMethod, aArg);
		// aMethod is a string - the method to call in framescript
		// aCallback is a function - optional - it will be triggered when aMethod is done calling

		var cbid = null;
		if (typeof(aMethod) == 'number') {
			// this is a response to a callack waiting in framescript
			cbid = aMethod;
			aMethod = null;
		} else {
			if (aCallback) {
				cbid = this.nextcbid++;
				this.callbackReceptacle[cbid] = aCallback;
			}
		}

		// return;
		aMessageManager.sendAsyncMessage(aChannelId, {
			method: aMethod,
			arg: aArg,
			cbid
		});
	};
	this.callbackReceptacle = {};

	Services.mm.addMessageListener(aChannelId, this.listener);
}
// start - CommAPI for bootstrap-framescript - bootstrap side - cross-file-link55565665464644
// start - CommAPI for bootstrap-content - bootstrap side - cross-file-link0048958576532536411
// message method - postMessage - content is in-process-content-windows, transferring works
// there is a framescript version of this, because framescript cant get aPort1 and aPort2 so it has to create its own
function contentComm(aContentWindow, aPort1, aPort2, onHandshakeComplete) {
	// onHandshakeComplete is triggered when handshake is complete
	// when a new contentWindow creates a contentComm on contentWindow side, it requests whatever it needs on init, so i dont offer a onBeforeInit. I do offer a onHandshakeComplete which is similar to onAfterInit, but not exactly the same
	// no unregister for this really, as no listeners setup, to unregister you just need to GC everything, so just break all references to it

	var handshakeComplete = false; // indicates this.postMessage will now work i think. it might work even before though as the messages might be saved till a listener is setup? i dont know i should ask
	var scope = gBootstrap;

	this.CallbackTransferReturn = function(aArg, aTransfers) {
		// aTransfers should be an array
		this.arg = aArg;
		this.xfer = aTransfers;
	};

	this.listener = function(e) {
		var payload = e.data;
		console.log('bootstrap contentComm - incoming, payload:', payload); //, 'e:', e);

		if (payload.method) {
			if (payload.method == 'contentComm_handshake_finalized') {
				handshakeComplete = false;
				if (onHandshakeComplete) {
					onHandshakeComplete(this);
				}
				return;
			}
			if (!(payload.method in scope)) { console.error('method of "' + payload.method + '" not in scope'); throw new Error('method of "' + payload.method + '" not in scope') } // dev line remove on prod
			var rez_bs_call_for_win = scope[payload.method](payload.arg, this);
			console.log('rez_bs_call_for_win:', rez_bs_call_for_win);
			if (payload.cbid) {
				if (rez_bs_call_for_win && rez_bs_call_for_win.constructor.name == 'Promise') {
					rez_bs_call_for_win.then(
						function(aVal) {
							console.log('Fullfilled - rez_bs_call_for_win - ', aVal);
							this.postMessage(payload.cbid, aVal);
						}.bind(this),
						genericReject.bind(null, 'rez_bs_call_for_win', 0)
					).catch(genericCatch.bind(null, 'rez_bs_call_for_win', 0));
				} else {
					console.log('calling postMessage for callback with rez_bs_call_for_win:', rez_bs_call_for_win, 'this:', this);
					this.postMessage(payload.cbid, rez_bs_call_for_win);
				}
			}
		} else if (!payload.method && payload.cbid) {
			// its a cbid
			this.callbackReceptacle[payload.cbid](payload.arg, this);
			delete this.callbackReceptacle[payload.cbid];
		} else {
			throw new Error('invalid combination');
		}
	}.bind(this);

	this.nextcbid = 1; //next callback id

	this.postMessage = function(aMethod, aArg, aTransfers, aCallback) {

		// aMethod is a string - the method to call in framescript
		// aCallback is a function - optional - it will be triggered when aMethod is done calling
		if (aArg && aArg.constructor == this.CallbackTransferReturn) {
			// aTransfers is undefined
			// i needed to create CallbackTransferReturn so that callbacks can transfer data back
			aTransfers = aArg.xfer;
			aArg = aArg.arg;
		}
		var cbid = null;
		if (typeof(aMethod) == 'number') {
			// this is a response to a callack waiting in framescript
			cbid = aMethod;
			aMethod = null;
		} else {
			if (aCallback) {
				cbid = this.nextcbid++;
				this.callbackReceptacle[cbid] = aCallback;
			}
		}

		// return;
		aPort1.postMessage({
			method: aMethod,
			arg: aArg,
			cbid
		}, aTransfers ? aTransfers : undefined);
	}

	aPort1.onmessage = this.listener;
	this.callbackReceptacle = {};

	aContentWindow.postMessage({
		topic: 'contentComm_handshake',
		port2: aPort2
	}, '*', [aPort2]);

}
// end - CommAPI for bootstrap-content - bootstrap side - cross-file-link0048958576532536411
// start - CommAPI for bootstrap-worker - bootstrap side - cross-file-link5323131347
// message method - postMessage
// on unregister, workers are terminated
var gWorkerComms = [];
function workerComm_unregAll() {
	var l = gWorkerComms.length;
	for (var i=0; i<l; i++) {
		gWorkerComms[i].unregister();
	}
}
function workerComm(aWorkerPath, onBeforeInit, onAfterInit, aWebWorker) {
	// limitations:
		// the first call is guranteed
		// devuser should never postMessage from worker with method name "triggerOnAfterInit" - this is reserved for programtic use
		// devuser should never postMessage from bootstrap with method name "init" - progmaticcaly this is automatically done in this.createWorker

	// worker is lazy loaded, it is not created until the first call. if you want instant instantiation, call this.createWorker() with no args
	// creates a ChromeWorker, unless aWebWorker is true

	// if onBeforeInit is set
		// if worker has `init` function
			// it is called by bootstrap, (progrmatically, i determine this by basing the first call to the worker)
	// if onBeforeInit is NOT set
		// if worker has `init` function
			// it is called by the worker before the first call to any method in the worker
	// onAfterInit is not called if `init` function does NOT exist in the worker. it is called by worker doing postMessage to bootstrap

	// onBeforeInit - args: this - it is a function, return a single var to send to init function in worker. can return CallbackTransferReturn if you want to transfer. it is run to build the data the worker should be inited with.
	// onAfterInit - args: aArg, this - a callback that happens after init is complete. aArg is return value of init from in worker. the first call to worker will happen after onAfterInit runs in bootstrap
	// these init features are offered because most times, workers need some data before starting off. and sometimes data is sent back to bootstrap like from init of MainWorker's
	// no featuere for prep term, as the prep term should be done in the `self.onclose = function(event) { ... }` of the worker
	gWorkerComms.push(this);

	var worker;
	var scope = gBootstrap;
	this.nextcbid = 1; //next callback id
	this.callbackReceptacle = {};
	this.CallbackTransferReturn = function(aArg, aTransfers) {
		// aTransfers should be an array
		this.arg = aArg;
		this.xfer = aTransfers;
	};
	this.createWorker = function(onAfterCreate) {
		// only triggered by postMessage when `var worker` has not yet been set
		worker = aWebWorker ? new Worker(aWorkerPath) : new ChromeWorker(aWorkerPath);
		worker.addEventListener('message', this.listener);

		if (onAfterInit) {
			var oldOnAfterInit = onAfterInit;
			onAfterInit = function(aArg, aComm) {
				oldOnAfterInit(aArg, aComm);
				if (onAfterCreate) {
					onAfterCreate(); // link39399999
				}
			}
		}

		var initArg;
		if (onBeforeInit) {
			initArg = onBeforeInit(this);
			this.postMessage('init', initArg); // i dont put onAfterCreate as a callback here, because i want to gurantee that the call of onAfterCreate happens after onAfterInit is triggered link39399999
		} else {
			// else, worker is responsible for calling init. worker will know because it keeps track in listener, what is the first postMessage, if it is not "init" then it will run init
			if (onAfterCreate) {
				onAfterCreate(); // as postMessage i the only one who calls this.createWorker(), onAfterCreate is the origianl postMessage intended by the devuser
			}
		}
	};
	this.postMessage = function(aMethod, aArg, aTransfers, aCallback) {
		// aMethod is a string - the method to call in framescript
		// aCallback is a function - optional - it will be triggered when aMethod is done calling

		if (!worker) {
			this.createWorker(this.postMessage.bind(this, aMethod, aArg, aTransfers, aCallback));
		} else {
			if (aArg && aArg.constructor == this.CallbackTransferReturn) {
				// aTransfers is undefined
				// i needed to create CallbackTransferReturn so that callbacks can transfer data back
				aTransfers = aArg.xfer;
				aArg = aArg.arg;
			}
			var cbid = null;
			if (typeof(aMethod) == 'number') {
				// this is a response to a callack waiting in framescript
				cbid = aMethod;
				aMethod = null;
			} else {
				if (aCallback) {
					cbid = this.nextcbid++;
					this.callbackReceptacle[cbid] = aCallback;
				}
			}

			worker.postMessage({
				method: aMethod,
				arg: aArg,
				cbid
			}, aTransfers ? aTransfers : undefined);
		}
	};
	this.unregister = function() {

		var l = gWorkerComms.length;
		for (var i=0; i<l; i++) {
			if (gWorkerComms[i] == this) {
				gWorkerComms.splice(i, 1);
				break;
			}
		}

		if (worker) { // because maybe it was setup, but never instantiated
			worker.terminate();
		}

	};
	this.listener = function(e) {
		var payload = e.data;
		console.log('bootstrap workerComm - incoming, payload:', payload); //, 'e:', e);

		if (payload.method) {
			if (payload.method == 'triggerOnAfterInit') {
				if (onAfterInit) {
					onAfterInit(payload.arg, this);
				}
				return;
			}
			if (!(payload.method in scope)) { console.error('method of "' + payload.method + '" not in scope'); throw new Error('method of "' + payload.method + '" not in scope') } // dev line remove on prod
			var rez_bs_call_for_worker = scope[payload.method](payload.arg, this);
			console.log('rez_bs_call_for_worker:', rez_bs_call_for_worker);
			if (payload.cbid) {
				if (rez_bs_call_for_worker && rez_bs_call_for_worker.constructor.name == 'Promise') {
					rez_bs_call_for_worker.then(
						function(aVal) {
							console.log('Fullfilled - rez_bs_call_for_worker - ', aVal);
							this.postMessage(payload.cbid, aVal);
						}.bind(this),
						genericReject.bind(null, 'rez_bs_call_for_worker', 0)
					).catch(genericCatch.bind(null, 'rez_bs_call_for_worker', 0));
				} else {
					console.log('calling postMessage for callback with rez_bs_call_for_worker:', rez_bs_call_for_worker, 'this:', this);
					this.postMessage(payload.cbid, rez_bs_call_for_worker);
				}
			}
		} else if (!payload.method && payload.cbid) {
			// its a cbid
			this.callbackReceptacle[payload.cbid](payload.arg, this);
			delete this.callbackReceptacle[payload.cbid];
		} else {
			console.error('bootstrap workerComm - invalid combination');
			throw new Error('bootstrap workerComm - invalid combination');
		}
	}.bind(this);
}
// end - CommAPI for bootstrap-worker - bootstrap side - cross-file-link5323131347
// end - CommAPI

// end - common helper functions
