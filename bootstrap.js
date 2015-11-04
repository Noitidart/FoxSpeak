// Imports
const {classes: Cc, interfaces: Ci, utils: Cu} = Components;
Cu.import('resource:///modules/CustomizableUI.jsm');
Cu.import('resource://gre/modules/devtools/Console.jsm');
Cu.import('resource://gre/modules/osfile.jsm');
Cu.import('resource://gre/modules/Promise.jsm');
Cu.import('resource://gre/modules/Services.jsm');
Cu.import('resource://gre/modules/XPCOMUtils.jsm');
// Cu.importGlobalProperties(['btoa']);

const PromiseWorker = Cu.import('resource://gre/modules/PromiseWorker.jsm', {}).BasePromiseWorker;


// Globals
const core = {
	addon: {
		name: 'FoxSpeak',
		id: 'FoxSpeak@jetpack',
		path: {
			name: 'foxspeak',
			content: 'chrome://foxspeak/content/',
			images: 'chrome://foxspeak/content/resources/images/',
			locale: 'chrome://foxspeak/locale/',
			modules: 'chrome://foxspeak/content/modules/',
			resources: 'chrome://foxspeak/content/resources/',
			styles: 'chrome://foxspeak/content/resources/styles/',
			workers: 'chrome://foxspeak/content/modules/workers/'
		},
		cache_key: Math.random() // set to version on release
	},
	os: {
		name: OS.Constants.Sys.Name.toLowerCase(),
		toolkit: Services.appinfo.widgetToolkit.toLowerCase(),
		xpcomabi: Services.appinfo.XPCOMABI
	}
};

const JETPACK_DIR_BASENAME = 'jetpack';
const OSPath_simpleStorage = OS.Path.join(OS.Constants.Path.profileDir, JETPACK_DIR_BASENAME, core.addon.id, 'simple-storage');
const myPrefBranch = 'extensions.' + core.addon.id + '.';

var bootstrap = this;
const NS_HTML = 'http://www.w3.org/1999/xhtml';
const CUI_CSSURI = Services.io.newURI(core.addon.path.styles + 'cui.css', null, null);

// Lazy Imports
const myServices = {};
XPCOMUtils.defineLazyGetter(myServices, 'hph', function () { return Cc['@mozilla.org/network/protocol;1?name=http'].getService(Ci.nsIHttpProtocolHandler); });
XPCOMUtils.defineLazyGetter(myServices, 'sb', function () { return Services.strings.createBundle(core.addon.path.locale + 'global.properties?' + Math.random()); /* Randomize URI to work around bug 719376 */ });

// START - Addon Functionalities
var myWidgetListener = {
	onWidgetDestroyed: function(aWidgetId) {
		if (aWidgetId != 'cui_foxspeak') {
			return
		}
		console.log('my widget destoryed');
		CustomizableUI.removeListener(myWidgetListener);
		
		var DOMWindows = Services.wm.getEnumerator('navigator:browser');
		while (DOMWindows.hasMoreElements()) {
			var DOMWindow = DOMWindows.getNext();
			var DOMWindowUtils = DOMWindow.QueryInterface(Ci.nsIInterfaceRequestor).getInterface(Ci.nsIDOMWindowUtils);
			DOMWindowUtils.removeSheet(CUI_CSSURI, DOMWindowUtils.AUTHOR_SHEET);
		}
		
		console.log('all my css styles removed from windows that had a cui, which are navigator:browser type windows');
	}
};
// END - Addon Functionalities

function install() {}

function uninstall(aData, aReason) {
	if (aReason == ADDON_UNINSTALL) {
		// delete prefs
	}
}

function startup(aData, aReason) {
	// core.addon.aData = aData;
	extendCore();
	
	// core.addon.path.jar = 'jar:' + OS.Path.toFileURI(aData.installPath.path).replace(aData.id + '.xpi', aData.id + '.xpi!/');
	core.addon.path.jar = 'jar:' + OS.Path.toFileURI(aData.installPath.path) + '!/';
	core.addon.path.file = aData.installPath.path;
	// jar:file:///C:/Users/Vayeate/AppData/Roaming/Mozilla/Firefox/Profiles/aksozfjt.Unnamed%20Profile%2010/extensions/FoxSpeak@jetpack.xpi!/
	// path to bootstrap.js would be: `jar:file:///C:/Users/Vayeate/AppData/Roaming/Mozilla/Firefox/Profiles/aksozfjt.Unnamed%20Profile%2010/extensions/FoxSpeak@jetpack.xpi!/bootstrap.js`
	
	console.log('core.addon.path.jar', core.addon.path.jar);
	
	var promise_getMainWorker = SIPWorker('MainWorker', core.addon.path.workers + 'MainWorker.js');
	promise_getMainWorker.then(
		function(aVal) {
			console.log('Fullfilled - promise_getMainWorker - ', aVal);
			// start - do stuff here - promise_getMainWorker
			
				// setup custom MainWorker messages - PromiseWorker style
				// Define a custom error prototype.
				function MainWorkerError(errObj) {
				  this.message = errObj.message;
				  this.name = errObj.name;
				}
				MainWorkerError.fromMsg = function(msgObj) {
				  return new MainWorkerError(msgObj);
				};

				// Register a constructor.
				MainWorker.ExceptionHandlers['MainWorkerError'] = MainWorkerError.fromMsg;
			
			// end - do stuff here - promise_getMainWorker
		},
		function(aReason) {
			var rejObj = {
				name: 'promise_getMainWorker',
				aReason: aReason
			};
			console.warn('Rejected - promise_getMainWorker - ', rejObj);
		}
	).catch(
		function(aCaught) {
			var rejObj = {
				name: 'promise_getMainWorker',
				aCaught: aCaught
			};
			console.error('Caught - promise_getMainWorker - ', rejObj);
		}
	);
	
	CustomizableUI.addListener(myWidgetListener);
	CustomizableUI.createWidget({
		id: 'cui_foxspeak',
		defaultArea: CustomizableUI.AREA_NAVBAR,
		label: myServices.sb.GetStringFromName('cui_foxspeak_lbl'),
		tooltiptext: myServices.sb.GetStringFromName('cui_foxspeak_tip'),
		onBeforeCreated: function(aDOMDocument) {
			console.error('onBeforeCreated triggered with aDOMDocument:', aDOMDocument);
			var DOMWindow = aDOMDocument.defaultView;
			var DOMWindowUtils = DOMWindow.QueryInterface(Ci.nsIInterfaceRequestor).getInterface(Ci.nsIDOMWindowUtils);
			console.log('CUI_CSSURI:', CUI_CSSURI);
			DOMWindowUtils.loadSheet(CUI_CSSURI, DOMWindowUtils.AUTHOR_SHEET);
			console.error('ok added my css');
		},
		onCommand: function(aEvent) {
			var DOMWindow = aEvent.target.ownerDocument.defaultView;
			
			DOMWindow.alert('hello from foxspeak!');
		}
	});
}

function shutdown(aData, aReason) {
	if (aReason == APP_SHUTDOWN) { return }
	
	CustomizableUI.destroyWidget('cui_foxspeak');
}

// start - common helper functions
function Deferred() {
	if (Promise && Promise.defer) {
		//need import of Promise.jsm for example: Cu.import('resource:/gree/modules/Promise.jsm');
		return Promise.defer();
	} else if (PromiseUtils && PromiseUtils.defer) {
		//need import of PromiseUtils.jsm for example: Cu.import('resource:/gree/modules/PromiseUtils.jsm');
		return PromiseUtils.defer();
	} else if (Promise) {
		try {
			/* A method to resolve the associated Promise with the value passed.
			 * If the promise is already settled it does nothing.
			 *
			 * @param {anything} value : This value is used to resolve the promise
			 * If the value is a Promise then the associated promise assumes the state
			 * of Promise passed as value.
			 */
			this.resolve = null;

			/* A method to reject the assocaited Promise with the value passed.
			 * If the promise is already settled it does nothing.
			 *
			 * @param {anything} reason: The reason for the rejection of the Promise.
			 * Generally its an Error object. If however a Promise is passed, then the Promise
			 * itself will be the reason for rejection no matter the state of the Promise.
			 */
			this.reject = null;

			/* A newly created Pomise object.
			 * Initially in pending state.
			 */
			this.promise = new Promise(function(resolve, reject) {
				this.resolve = resolve;
				this.reject = reject;
			}.bind(this));
			Object.freeze(this);
		} catch (ex) {
			console.error('Promise not available!', ex);
			throw new Error('Promise not available!');
		}
	} else {
		throw new Error('Promise not available!');
	}
}

function SIPWorker(workerScopeName, aPath, aCore=core) {
	// "Start and Initialize PromiseWorker"
	// returns promise
		// resolve value: jsBool true
	// aCore is what you want aCore to be populated with
	// aPath is something like `core.addon.path.content + 'modules/workers/blah-blah.js'`
	
	// :todo: add support and detection for regular ChromeWorker // maybe? cuz if i do then ill need to do ChromeWorker with callback
	
	var deferredMain_SIPWorker = new Deferred();

	if (!(workerScopeName in bootstrap)) {
		bootstrap[workerScopeName] = new PromiseWorker(aPath);
		
		if ('addon' in aCore && 'aData' in aCore.addon) {
			delete aCore.addon.aData; // we delete this because it has nsIFile and other crap it, but maybe in future if I need this I can try JSON.stringify'ing it
		}
		
		var promise_initWorker = bootstrap[workerScopeName].post('init', [aCore]);
		promise_initWorker.then(
			function(aVal) {
				console.log('Fullfilled - promise_initWorker - ', aVal);
				// start - do stuff here - promise_initWorker
				deferredMain_SIPWorker.resolve(true);
				// end - do stuff here - promise_initWorker
			},
			function(aReason) {
				var rejObj = {name:'promise_initWorker', aReason:aReason};
				console.warn('Rejected - promise_initWorker - ', rejObj);
				deferredMain_SIPWorker.reject(rejObj);
			}
		).catch(
			function(aCaught) {
				var rejObj = {name:'promise_initWorker', aCaught:aCaught};
				console.error('Caught - promise_initWorker - ', rejObj);
				deferredMain_SIPWorker.reject(rejObj);
			}
		);
		
	} else {
		deferredMain_SIPWorker.reject('Something is loaded into bootstrap[workerScopeName] already');
	}
	
	return deferredMain_SIPWorker.promise;
	
}

function aReasonMax(aReason) {
	var deepestReason = aReason;
	while (deepestReason.hasOwnProperty('aReason') || deepestReason.hasOwnProperty()) {
		if (deepestReason.hasOwnProperty('aReason')) {
			deepestReason = deepestReason.aReason;
		} else if (deepestReason.hasOwnProperty('aCaught')) {
			deepestReason = deepestReason.aCaught;
		}
	}
	return deepestReason;
}

function extendCore() {
	// adds some properties i use to core based on the current operating system, it needs a switch, thats why i couldnt put it into the core obj at top
	switch (core.os.name) {
		case 'winnt':
		case 'winmo':
		case 'wince':
			core.os.version = parseFloat(Services.sysinfo.getProperty('version'));
			// http://en.wikipedia.org/wiki/List_of_Microsoft_Windows_versions
			if (core.os.version == 6.0) {
				core.os.version_name = 'vista';
			}
			if (core.os.version >= 6.1) {
				core.os.version_name = '7+';
			}
			if (core.os.version == 5.1 || core.os.version == 5.2) { // 5.2 is 64bit xp
				core.os.version_name = 'xp';
			}
			break;
			
		case 'darwin':
			var userAgent = myServices.hph.userAgent;

			var version_osx = userAgent.match(/Mac OS X 10\.([\d\.]+)/);

			
			if (!version_osx) {
				throw new Error('Could not identify Mac OS X version.');
			} else {
				var version_osx_str = version_osx[1];
				var ints_split = version_osx[1].split('.');
				if (ints_split.length == 1) {
					core.os.version = parseInt(ints_split[0]);
				} else if (ints_split.length >= 2) {
					core.os.version = ints_split[0] + '.' + ints_split[1];
					if (ints_split.length > 2) {
						core.os.version += ints_split.slice(2).join('');
					}
					core.os.version = parseFloat(core.os.version);
				}
				// this makes it so that 10.10.0 becomes 10.100
				// 10.10.1 => 10.101
				// so can compare numerically, as 10.100 is less then 10.101
				
				//core.os.version = 6.9; // note: debug: temporarily forcing mac to be 10.6 so we can test kqueue
			}
			break;
		default:
			// nothing special
	}
	

}
// end - common helper functions