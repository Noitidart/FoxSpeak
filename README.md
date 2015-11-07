# FoxSpeak
ff-addon: An accessibility add-on which allows voice-to-text and text-to-voice.

The API being used is the pocketsphinx API by CMU. (http://cmusphinx.sourceforge.net/doc/pocketsphinx/)
Currently work is being done on the linux addon.

To do :
1) A .dll using the functions defined in pocketsphinx for windows,
2) A .dylib for iOS, so that @noitidart can start writing ctypes for both.
3) Port this on android.
4) Decide interface, etc. (of course after the backend work is done).
