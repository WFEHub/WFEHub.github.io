/*
    PRIZM address class, extended version (with error guessing).

    Version: 1.0, license: Public Domain, coder: modified PRIZM GROUP.
*/

function PrizmAddress() {
	var codeword = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
	var syndrome = [0, 0, 0, 0, 0];

	var gexp = [1, 2, 4, 8, 16, 5, 10, 20, 13, 26, 17, 7, 14, 28, 29, 31, 27, 19, 3, 6, 12, 24, 21, 15, 30, 25, 23, 11, 22, 9, 18, 1];
	var glog = [0, 0, 1, 18, 2, 5, 19, 11, 3, 29, 6, 27, 20, 8, 12, 23, 4, 10, 30, 17, 7, 22, 28, 26, 21, 25, 9, 16, 13, 14, 24, 15];

	var cwmap = [3, 2, 1, 0, 7, 6, 5, 4, 13, 14, 15, 16, 12, 8, 9, 10, 11];

	var alphabet = 'PRZM23456789ABCDEFGHJKLNQSTUVWXY';

	this.guess = [];

	function ginv(a) {
		return gexp[31 - glog[a]];
	}

	function gmult(a, b) {
		if (a == 0 || b == 0) return 0;

		var idx = (glog[a] + glog[b]) % 31;

		return gexp[idx];
	} //__________________________

	function calc_discrepancy(lambda, r) {
		var discr = 0;

		for (var i = 0; i < r; i++) {
			discr ^= gmult(lambda[i], syndrome[r - i]);
		}

		return discr;
	} //__________________________

	function find_errors(lambda) {
		var errloc = [];

		for (var i = 1; i <= 31; i++) {
			var sum = 0;

			for (var j = 0; j < 5; j++) {
				sum ^= gmult(gexp[(j * i) % 31], lambda[j]);
			}

			if (sum == 0) {
				var pos = 31 - i;
				if (pos > 12 && pos < 27) return [];

				errloc[errloc.length] = pos;
			}
		}

		return errloc;
	} //__________________________

	function guess_errors() {
		var el = 0,
			b = [0, 0, 0, 0, 0],
			t = [];

		var deg_lambda = 0,
			lambda = [1, 0, 0, 0, 0]; // error+erasure locator poly

		// Berlekamp-Massey algorithm to determine error+erasure locator polynomial

		for (var r = 0; r < 4; r++) {
			var discr = calc_discrepancy(lambda, r + 1); // Compute discrepancy at the r-th step in poly-form

			if (discr != 0) {
				deg_lambda = 0;

				for (var i = 0; i < 5; i++) {
					t[i] = lambda[i] ^ gmult(discr, b[i]);

					if (t[i]) deg_lambda = i;
				}

				if (2 * el <= r) {
					el = r + 1 - el;

					for (i = 0; i < 5; i++) {
						b[i] = gmult(lambda[i], ginv(discr));
					}
				}

				lambda = t.slice(); // copy
			}

			b.unshift(0); // shift => mul by x
		}

		// Find roots of the locator polynomial.

		var errloc = find_errors(lambda);

		var errors = errloc.length;

		if (errors < 1 || errors > 2) return false;

		if (deg_lambda != errors) return false; // deg(lambda) unequal to number of roots => uncorrectable error

		// Compute err+eras evaluator poly omega(x) = s(x)*lambda(x) (modulo x**(4)). Also find deg(omega).

		var omega = [0, 0, 0, 0, 0];

		for (var i = 0; i < 4; i++) {
			var t = 0;

			for (var j = 0; j < i; j++) {
				t ^= gmult(syndrome[i + 1 - j], lambda[j]);
			}

			omega[i] = t;
		}

		// Compute error values in poly-form.

		for (r = 0; r < errors; r++) {
			var t = 0;
			var pos = errloc[r];
			var root = 31 - pos;

			for (i = 0; i < 4; i++) // evaluate Omega at alpha^(-i)
			{
				t ^= gmult(omega[i], gexp[(root * i) % 31]);
			}

			if (t) // evaluate Lambda' (derivative) at alpha^(-i); all odd powers disappear
			{
				var denom = gmult(lambda[1], 1) ^ gmult(lambda[3], gexp[(root * 2) % 31]);

				if (denom == 0) return false;

				if (pos > 12) pos -= 14;

				codeword[pos] ^= gmult(t, ginv(denom));
			}
		}

		return true;
	} //__________________________

	function encode() {
		var p = [0, 0, 0, 0];

		for (var i = 12; i >= 0; i--) {
			var fb = codeword[i] ^ p[3];

			p[3] = p[2] ^ gmult(30, fb);
			p[2] = p[1] ^ gmult(6, fb);
			p[1] = p[0] ^ gmult(9, fb);
			p[0] = gmult(17, fb);
		}

		codeword[13] = p[0];
		codeword[14] = p[1];
		codeword[15] = p[2];
		codeword[16] = p[3];
	} //__________________________

	function reset() {
		for (var i = 0; i < 17; i++) codeword[i] = 1;
	} //__________________________

	function set_codeword(cw, len, skip) {
		if (typeof len === 'undefined') len = 17;
		if (typeof skip === 'undefined') skip = -1;

		for (var i = 0, j = 0; i < len; i++) {
			if (i != skip) codeword[cwmap[j++]] = cw[i];
		}
	} //__________________________

	this.add_guess = function() {
		var s = this.toString(),
			len = this.guess.length;

		if (len > 2) return;

		for (var i = 0; i < len; i++) {
			if (this.guess[i] == s) return;
		}

		this.guess[len] = s;
	} //__________________________

	this.ok = function() {
		var sum = 0;

		for (var i = 1; i < 5; i++) {
			for (var j = 0, t = 0; j < 31; j++) {
				if (j > 12 && j < 27) continue;

				var pos = j;
				if (j > 26) pos -= 14;

				t ^= gmult(codeword[pos], gexp[(i * j) % 31]);
			}

			sum |= t;
			syndrome[i] = t;
		}

		return (sum == 0);
	} //__________________________

	function from_acc(acc) {
		var inp = [],
			out = [],
			pos = 0,
			len = acc.length;

		if (len == 20 && acc.charAt(0) != '1') return false;

		for (var i = 0; i < len; i++) {
			inp[i] = acc.charCodeAt(i) - '0'.charCodeAt(0);
		}

		do // base 10 to base 32 conversion
		{
			var divide = 0,
				newlen = 0;

			for (i = 0; i < len; i++) {
				divide = divide * 10 + inp[i];

				if (divide >= 32) {
					inp[newlen++] = divide >> 5;
					divide &= 31;
				} else if (newlen > 0) {
					inp[newlen++] = 0;
				}
			}

			len = newlen;
			out[pos++] = divide;
		}
		while (newlen);

		for (i = 0; i < 13; i++) // copy to codeword in reverse, pad with 0's
		{
			codeword[i] = (--pos >= 0 ? out[i] : 0);
		}

		encode();

		return true;
	} //__________________________

	this.toString = function() {
		var out = 'PRIZM-';

		for (var i = 0; i < 17; i++) {
			out += alphabet[codeword[cwmap[i]]];

			if ((i & 3) == 3 && i < 13) out += '-';
		}

		return out;
	} //__________________________

	this.account_id = function() {
		var out = '',
			inp = [],
			len = 13;

		for (var i = 0; i < 13; i++) {
			inp[i] = codeword[12 - i];
		}

		do // base 32 to base 10 conversion
		{
			var divide = 0,
				newlen = 0;

			for (i = 0; i < len; i++) {
				divide = divide * 32 + inp[i];

				if (divide >= 10) {
					inp[newlen++] = Math.floor(divide / 10);
					divide %= 10;
				} else if (newlen > 0) {
					inp[newlen++] = 0;
				}
			}

			len = newlen;
			out += String.fromCharCode(divide + '0'.charCodeAt(0));
		}
		while (newlen);

		return out.split("").reverse().join("");
	} //__________________________

	this.set = function(adr, allow_accounts) {
		if (typeof allow_accounts === 'undefined') allow_accounts = true;

		var len = 0;
		this.guess = [];
		reset();

		adr = String(adr);

		adr = adr.replace(/(^\s+)|(\s+$)/g, '').toUpperCase();

		if (adr.indexOf('PRIZM-') == 0) adr = adr.substr(6);

		if (adr.match(/^\d{1,20}$/g)) // account id
		{
			if (allow_accounts) return from_acc(adr);
		} else // address
		{
			var clean = [];

			for (var i = 0; i < adr.length; i++) {
				var pos = alphabet.indexOf(adr[i]);

				if (pos >= 0) {
					clean[len++] = pos;
					if (len > 18) return false;
				}
			}
		}

		if (len == 16) // guess deletion
		{
			for (var i = 16; i >= 0; i--) {
				for (var j = 0; j < 32; j++) {
					clean[i] = j;

					set_codeword(clean);

					if (this.ok()) this.add_guess();
				}

				if (i > 0) {
					var t = clean[i - 1];
					clean[i - 1] = clean[i];
					clean[i] = t;
				}
			}
		}

		if (len == 18) // guess insertion
		{
			for (var i = 0; i < 18; i++) {
				set_codeword(clean, 18, i);

				if (this.ok()) this.add_guess();
			}
		}

		if (len == 17) {
			set_codeword(clean);

			if (this.ok()) return true;

			if (guess_errors() && this.ok()) this.add_guess();
		}

		reset();

		return false;
	}

	this.format_guess = function(s, org) {
		var d = '',
			list = [];

		s = s.toUpperCase();
		org = org.toUpperCase();

		for (var i = 0; i < s.length;) {
			var m = 0;

			for (var j = 1; j < s.length; j++) {
				var pos = org.indexOf(s.substr(i, j));

				if (pos != -1) {
					if (Math.abs(pos - i) < 3) m = j;
				} else break;
			}

			if (m) {
				list[list.length] = {
					's': i,
					'e': i + m
				};
				i += m;
			} else i++;
		}

		if (list.length == 0) return s;

		for (var i = 0, j = 0; i < s.length; i++) {
			if (i >= list[j].e) {
				var start;

				while (j < list.length - 1) {
					start = list[j++].s;

					if (i < list[j].e || list[j].s >= start) break;
				}
			}

			if (i >= list[j].s && i < list[j].e) {
				d += s.charAt(i);
			} else {
				d += '<b style="color:red">' + s.charAt(i) + '</b>';
			}
		}

		return d;
	};
}

/* global converters, curve25519, CryptoJS */

function signBytes(message, secretPhrase) {
    if (!secretPhrase) {
        throw {"message": $.t("error_encryption_passphrase_required"), "errorCode": 1};
    }
    var messageBytes = converters.hexStringToByteArray(message);
    var secretPhraseBytes = converters.hexStringToByteArray(converters.stringToHexString(secretPhrase));

    var digest = simpleHash(secretPhraseBytes);
    var s = curve25519.keygen(digest).s;
    var m = simpleHash(messageBytes);
    var x = simpleHash(m, s);
    var y = curve25519.keygen(x).p;
    var h = simpleHash(m, y);
    var v = curve25519.sign(h, x, s);
//		return converters.byteArrayToHexString(v.concat(h));
    var signature = converters.byteArrayToHexString(v.concat(h));

    var payload = message.substr(0, 192) + signature + message.substr(320);

    return payload;
}

function simpleHash(b1, b2) {
    var sha256 = CryptoJS.algo.SHA256.create();
    sha256.update(converters.byteArrayToWordArray(b1));
    if (b2) {
        sha256.update(converters.byteArrayToWordArray(b2));
    }
    var hash = sha256.finalize();
    return converters.wordArrayToByteArrayImpl(hash, false);
}

function gAIbyteArrayToBigInteger(byteArray) {
    var value = new BigInteger("0", 10);
    var temp1, temp2;
    for (var i = byteArray.length - 1; i >= 0; i--) {
        temp1 = value.multiply(new BigInteger("256", 10));
        temp2 = temp1.add(new BigInteger(byteArray[i].toString(10), 10));
        value = temp2;
    }

    return value;
}

function getPublicKeyPrizm(secretPhrase) {
    var secretPhraseBytes = converters.hexStringToByteArray(converters.stringToHexString(secretPhrase));
    var digest = simpleHash(secretPhraseBytes);
    return converters.byteArrayToHexString(curve25519.keygen(digest).p);
}

function getPrivateKey(secretPhrase) {
    var bytes = simpleHash(converters.stringToByteArray(secretPhrase));
    return converters.shortArrayToHexString(curve25519_clamp(converters.byteArrayToShortArray(bytes)));
    
}

function getAccountId(publicKey) {
    try {
        var hex = converters.hexStringToByteArray(publicKey);
        var account = simpleHash(hex);
        account = converters.byteArrayToHexString(account);
        var slice = (converters.hexStringToByteArray(account)).slice(0, 8);
        var accountId = gAIbyteArrayToBigInteger(slice).toString();
        return accountId;
    } catch (err) {
        return null;
    }

}

function getRSaddressPrizm(accountId) {
    var address = new PrizmAddress();
    if (address.set(accountId)) {
        return address.toString();
    } else {
        return "no address";
    }
}

function getIDByRSaddressPrizm(accountAddress) {
    var address = new PrizmAddress();
    if (address.set(accountAddress)) {
        return address.account_id();
    } else {
        return "no address";
    }
}

function getSharedSecret(key1, key2) {
    return converters.shortArrayToByteArray(curve25519_(converters.byteArrayToShortArray(key1), converters.byteArrayToShortArray(key2), null));
}

function aesDecrypt(ivCiphertext, options) {
    if (ivCiphertext.length < 16 || ivCiphertext.length % 16 !== 0) {
        throw {
            name: "invalid ciphertext"
        };
    }

    var iv = converters.byteArrayToWordArray(ivCiphertext.slice(0, 16));
    var ciphertext = converters.byteArrayToWordArray(ivCiphertext.slice(16));
    var sharedKey;
    if (!options.sharedKey) {
        sharedKey = getSharedSecret(options.privateKey, options.publicKey);
    } else {
        sharedKey = options.sharedKey.slice(0); //clone
    }

    var key;
    if (options.nonce) {
        for (var i = 0; i < 32; i++) {
            sharedKey[i] ^= options.nonce[i];
        }
        key = CryptoJS.SHA256(converters.byteArrayToWordArray(sharedKey));
    } else {
        key = converters.byteArrayToWordArray(sharedKey);
    }

    var encrypted = CryptoJS.lib.CipherParams.create({
        ciphertext: ciphertext,
        iv: iv,
        key: key
    });

    var decrypted = CryptoJS.AES.decrypt(encrypted, key, {
        iv: iv
    });

    return converters.wordArrayToByteArray(decrypted);
}

function aesEncrypt(plaintext, options) {
    if (!window.crypto && !window.msCrypto) {
        throw {
            "errorCode": -1,
            "message": "error_encryption_browser_support"
        };
    }

    // CryptoJS likes WordArray parameters
    var text = converters.byteArrayToWordArray(plaintext);
    var sharedKey;
    if (!options.sharedKey) {
        sharedKey = getSharedSecret(options.privateKey, options.publicKey);
    } else {
        sharedKey = options.sharedKey.slice(0); //clone
    }

    for (var i = 0; i < 32; i++) {
        sharedKey[i] ^= options.nonce[i];
    }

    var key = CryptoJS.SHA256(converters.byteArrayToWordArray(sharedKey));

    var tmp = new Uint8Array(16);

    if (window.crypto) {
        window.crypto.getRandomValues(tmp);
    } else {
        window.msCrypto.getRandomValues(tmp);
    }

    var iv = converters.byteArrayToWordArray(tmp);
    var encrypted = CryptoJS.AES.encrypt(text, key, {
        iv: iv
    });

    var ivOut = converters.wordArrayToByteArray(encrypted.iv);

    var ciphertextOut = converters.wordArrayToByteArray(encrypted.ciphertext);

    return ivOut.concat(ciphertextOut);
}

function decryptData(data, options) {
    if (!options.sharedKey) {
        options.sharedKey = getSharedSecret(options.privateKey, options.publicKey);
    }

    var compressedPlaintext = aesDecrypt(data, options);
    var binData = new Uint8Array(compressedPlaintext);
    return converters.byteArrayToString(pako.inflate(binData));
}

function encryptData(plaintext, options) {
    if (!window.crypto && !window.msCrypto) {
        throw {
            "errorCode": -1,
            "message": "error_encryption_browser_support"
        };
    }

    if (!options.sharedKey) {
        options.sharedKey = getSharedSecret(options.privateKey, options.publicKey);
    }

    var compressedPlaintext = pako.gzip(new Uint8Array(plaintext));
    options.nonce = new Uint8Array(32);
    if (window.crypto) {
//noinspection JSUnresolvedFunction
        window.crypto.getRandomValues(options.nonce);
    } else {
//noinspection JSUnresolvedFunction
        window.msCrypto.getRandomValues(options.nonce);
    }

    var data = aesEncrypt(compressedPlaintext, options);
    return {
        "nonce": options.nonce,
        "data": data
    };
}

function decryptNote(message, options, secretPhrase) {
    try {
        if (!options.sharedKey) {
            if (!options.privateKey) {
                options.privateKey = converters.hexStringToByteArray(getPrivateKey(secretPhrase));
            }

            if (!options.publicKey) {
                if (!options.account) {
                    throw {
                        "message": $.t("error_account_id_not_specified"),
                        "errorCode": 2
                    };
                }

                options.publicKey = converters.hexStringToByteArray(getPublicKeyPrizm(secretPhrase));
            }
        }

        options.nonce = converters.hexStringToByteArray(options.nonce);

        return decryptData(converters.hexStringToByteArray(message), options);
    } catch (err) {
        if (err.errorCode && err.errorCode < 3) {
            throw err;
        } else {
            throw {
                "message": "error_message_decryption",
                "errorCode": 3
            };
        }
    } 
}

function encryptNote(message, options, secretPhrase) {
    try {
        if (!options.sharedKey) {
            if (!options.privateKey) {
                if (!secretPhrase) {
                        throw {
                            "message": "error_encryption_passphrase_required",
                            "errorCode": 1
                        };
                }

                options.privateKey = converters.hexStringToByteArray(getPrivateKey(secretPhrase));
            }

            if (!options.publicKey) {
                if (!options.account) {
                    throw {
                        "message": "error_account_id_not_specified",
                        "errorCode": 2
                    };
                }

                try {
                    options.publicKey = converters.hexStringToByteArray(getPublicKeyPrizm(options.account, true));
                } catch (err) {
                    var pzmAddress = new PrizmAddress();

                    if (!pzmAddress.set(options.account)) {
                        throw {
                            "message": "error_invalid_account_id",
                            "errorCode": 3
                        };
                    } else {
                        throw {
                            "message": "error_public_key_not_specified",
                            "errorCode": 4
                        };
                    }
                }
            } else if (typeof options.publicKey == "string") {
                options.publicKey = converters.hexStringToByteArray(options.publicKey);
            }
        }

        var encrypted = encryptData(converters.stringToByteArray(message), options);

        return {
            "message": converters.byteArrayToHexString(encrypted.data),
            "nonce": converters.byteArrayToHexString(encrypted.nonce)
        };
    } catch (err) {
        return err;
    }
}
;

function decMessage(encrypted, nonce, publicKey, password) {
    try {
        var data = decryptNote(encrypted, {
            "nonce": nonce,
            "publicKey": converters.hexStringToByteArray(publicKey)
        }, password);
        return data;
    } catch (err) {
        return "error";
    }
}

function encMessage(message, account, publicKey, password) {
    try {
        var data = encryptNote(message, {
            "account": account,
            "publicKey": converters.hexStringToByteArray(publicKey)
        }, password);
        return data["message"] + ":" + data["nonce"];
    } catch (err) {
        return "error";
    }
}

function qrCodemaster(id, message) {
    document.getElementById(id).innerHTML = "";
    var qrcode = new QRCode(document.getElementById(id), {
        width: 600,
        height: 600,
        colorLight: "#f5f5f5"
    });
    qrcode.makeCode(message);
}

/******************************************************************************
 * Copyright © 2013-2016 The PZM Developers.                             *
 *                                                                            *
 * See the AUTHORS.txt, DEVELOPER-AGREEMENT.txt and LICENSE.txt files at      *
 * the top-level directory of this distribution for the individual copyright  *
 * holder information and the developer policies on copyright and licensing.  *
 *                                                                            *
 * Unless otherwise agreed in a custom licensing agreement, no part of the    *
 * PZM software, including this file, may be copied, modified, propagated,    *
 * or distributed except according to the terms contained in the LICENSE.txt  *
 * file.                                                                      *
 *                                                                            *
 * Removal or modification of this copyright notice is prohibited.            *
 *                                                                            *
 ******************************************************************************/

var converters = function() {
	var charToNibble = {};
	var nibbleToChar = [];
	var i;
	for (i = 0; i <= 9; ++i) {
		var character = i.toString();
		charToNibble[character] = i;
		nibbleToChar.push(character);
	}

	for (i = 10; i <= 15; ++i) {
		var lowerChar = String.fromCharCode('a'.charCodeAt(0) + i - 10);
		var upperChar = String.fromCharCode('A'.charCodeAt(0) + i - 10);

		charToNibble[lowerChar] = i;
		charToNibble[upperChar] = i;
		nibbleToChar.push(lowerChar);
	}

	return {
		byteArrayToHexString: function(bytes) {
			var str = '';
			for (var i = 0; i < bytes.length; ++i) {
				if (bytes[i] < 0) {
					bytes[i] += 256;
				}
				str += nibbleToChar[bytes[i] >> 4] + nibbleToChar[bytes[i] & 0x0F];
			}

			return str;
		},
		stringToByteArray: function(str) {
			str = unescape(encodeURIComponent(str)); //temporary

			var bytes = new Array(str.length);
			for (var i = 0; i < str.length; ++i)
				bytes[i] = str.charCodeAt(i);

			return bytes;
		},
		hexStringToByteArray: function(str) {
			var bytes = [];
			var i = 0;
			if (0 !== str.length % 2) {
				bytes.push(charToNibble[str.charAt(0)]);
				++i;
			}

			for (; i < str.length - 1; i += 2)
				bytes.push((charToNibble[str.charAt(i)] << 4) + charToNibble[str.charAt(i + 1)]);

			return bytes;
		},
		stringToHexString: function(str) {
			return this.byteArrayToHexString(this.stringToByteArray(str));
		},
		hexStringToString: function(hex) {
			return this.byteArrayToString(this.hexStringToByteArray(hex));
		},
		checkBytesToIntInput: function(bytes, numBytes, opt_startIndex) {
			var startIndex = opt_startIndex || 0;
			if (startIndex < 0) {
				throw new Error('Start index should not be negative');
			}

			if (bytes.length < startIndex + numBytes) {
				throw new Error('Need at least ' + (numBytes) + ' bytes to convert to an integer');
			}
			return startIndex;
		},
		byteArrayToSignedShort: function(bytes, opt_startIndex) {
			var index = this.checkBytesToIntInput(bytes, 2, opt_startIndex);
			var value = bytes[index];
			value += bytes[index + 1] << 8;
			return value;
		},
		byteArrayToSignedInt32: function(bytes, opt_startIndex) {
			var index = this.checkBytesToIntInput(bytes, 4, opt_startIndex);
			value = bytes[index];
			value += bytes[index + 1] << 8;
			value += bytes[index + 2] << 16;
			value += bytes[index + 3] << 24;
			return value;
		},
		byteArrayToBigInteger: function(bytes, opt_startIndex) {
			var index = this.checkBytesToIntInput(bytes, 8, opt_startIndex);

			var value = new BigInteger("0", 10);

			var temp1, temp2;

			for (var i = 7; i >= 0; i--) {
				temp1 = value.multiply(new BigInteger("256", 10));
				temp2 = temp1.add(new BigInteger(bytes[opt_startIndex + i].toString(10), 10));
				value = temp2;
			}

			return value;
		},
		// create a wordArray that is Big-Endian
		byteArrayToWordArray: function(byteArray) {
			var i = 0,
				offset = 0,
				word = 0,
				len = byteArray.length;
			var words = new Uint32Array(((len / 4) | 0) + (len % 4 == 0 ? 0 : 1));

			while (i < (len - (len % 4))) {
				words[offset++] = (byteArray[i++] << 24) | (byteArray[i++] << 16) | (byteArray[i++] << 8) | (byteArray[i++]);
			}
			if (len % 4 != 0) {
				word = byteArray[i++] << 24;
				if (len % 4 > 1) {
					word = word | byteArray[i++] << 16;
				}
				if (len % 4 > 2) {
					word = word | byteArray[i++] << 8;
				}
				words[offset] = word;
			}
			var wordArray = new Object();
			wordArray.sigBytes = len;
			wordArray.words = words;

			return wordArray;
		},
		// assumes wordArray is Big-Endian
		wordArrayToByteArray: function(wordArray) {
			return converters.wordArrayToByteArrayImpl(wordArray, true);
		},
		wordArrayToByteArrayImpl: function(wordArray, isFirstByteHasSign) {
			var len = wordArray.words.length;
			if (len == 0) {
				return new Array(0);
			}
			var byteArray = new Array(wordArray.sigBytes);
			var offset = 0,
				word, i;
			for (i = 0; i < len - 1; i++) {
				word = wordArray.words[i];
				byteArray[offset++] = isFirstByteHasSign ? word >> 24 : (word >> 24) & 0xff;
				byteArray[offset++] = (word >> 16) & 0xff;
				byteArray[offset++] = (word >> 8) & 0xff;
				byteArray[offset++] = word & 0xff;
			}
			word = wordArray.words[len - 1];
			byteArray[offset++] = isFirstByteHasSign ? word >> 24 : (word >> 24) & 0xff;
			if (wordArray.sigBytes % 4 == 0) {
				byteArray[offset++] = (word >> 16) & 0xff;
				byteArray[offset++] = (word >> 8) & 0xff;
				byteArray[offset++] = word & 0xff;
			}
			if (wordArray.sigBytes % 4 > 1) {
				byteArray[offset++] = (word >> 16) & 0xff;
			}
			if (wordArray.sigBytes % 4 > 2) {
				byteArray[offset++] = (word >> 8) & 0xff;
			}
			return byteArray;
		},
		byteArrayToString: function(bytes, opt_startIndex, length) {
			if (length == 0) {
				return "";
			}

			if (opt_startIndex && length) {
				var index = this.checkBytesToIntInput(bytes, parseInt(length, 10), parseInt(opt_startIndex, 10));

				bytes = bytes.slice(opt_startIndex, opt_startIndex + length);
			}

			return decodeURIComponent(escape(String.fromCharCode.apply(null, bytes)));
		},
		byteArrayToShortArray: function(byteArray) {
			var shortArray = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
			var i;
			for (i = 0; i < 16; i++) {
				shortArray[i] = byteArray[i * 2] | byteArray[i * 2 + 1] << 8;
			}
			return shortArray;
		},
		shortArrayToByteArray: function(shortArray) {
			var byteArray = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
			var i;
			for (i = 0; i < 16; i++) {
				byteArray[2 * i] = shortArray[i] & 0xff;
				byteArray[2 * i + 1] = shortArray[i] >> 8;
			}

			return byteArray;
		},
		shortArrayToHexString: function(ary) {
			var res = "";
			for (var i = 0; i < ary.length; i++) {
				res += nibbleToChar[(ary[i] >> 4) & 0x0f] + nibbleToChar[ary[i] & 0x0f] + nibbleToChar[(ary[i] >> 12) & 0x0f] + nibbleToChar[(ary[i] >> 8) & 0x0f];
			}
			return res;
		},
		/**
		 * Produces an array of the specified number of bytes to represent the integer
		 * value. Default output encodes ints in little endian format. Handles signed
		 * as well as unsigned integers. Due to limitations in JavaScript's number
		 * format, x cannot be a true 64 bit integer (8 bytes).
		 */
		intToBytes_: function(x, numBytes, unsignedMax, opt_bigEndian) {
			var signedMax = Math.floor(unsignedMax / 2);
			var negativeMax = (signedMax + 1) * -1;
			if (x != Math.floor(x) || x < negativeMax || x > unsignedMax) {
				throw new Error(
					x + ' is not a ' + (numBytes * 8) + ' bit integer');
			}
			var bytes = [];
			var current;
			// Number type 0 is in the positive int range, 1 is larger than signed int,
			// and 2 is negative int.
			var numberType = x >= 0 && x <= signedMax ? 0 :
				x > signedMax && x <= unsignedMax ? 1 : 2;
			if (numberType == 2) {
				x = (x * -1) - 1;
			}
			for (var i = 0; i < numBytes; i++) {
				if (numberType == 2) {
					current = 255 - (x % 256);
				} else {
					current = x % 256;
				}

				if (opt_bigEndian) {
					bytes.unshift(current);
				} else {
					bytes.push(current);
				}

				if (numberType == 1) {
					x = Math.floor(x / 256);
				} else {
					x = x >> 8;
				}
			}
			return bytes;

		},
		int32ToBytes: function(x, opt_bigEndian) {
			return converters.intToBytes_(x, 4, 4294967295, opt_bigEndian);
		},
        /**
         * Based on https://groups.google.com/d/msg/crypto-js/TOb92tcJlU0/Eq7VZ5tpi-QJ
         * Converts a word array to a Uint8Array.
         * @param {WordArray} wordArray The word array.
         * @return {Uint8Array} The Uint8Array.
         */
        wordArrayToByteArrayEx: function (wordArray) {
            // Shortcuts
            var words = wordArray.words;
            var sigBytes = wordArray.sigBytes;

            // Convert
            var u8 = new Uint8Array(sigBytes);
            for (var i = 0; i < sigBytes; i++) {
                var byte = (words[i >>> 2] >>> (24 - (i % 4) * 8)) & 0xff;
                u8[i]=byte;
            }

            return u8;
        },
        /**
         * Converts a Uint8Array to a word array.
         * @param {string} u8Str The Uint8Array.
         * @return {WordArray} The word array.
         */
        byteArrayToWordArrayEx: function (u8arr) {
            // Shortcut
            var len = u8arr.length;

            // Convert
            var words = [];
            for (var i = 0; i < len; i++) {
                words[i >>> 2] |= (u8arr[i] & 0xff) << (24 - (i % 4) * 8);
            }

            return CryptoJS.lib.WordArray.create(words, len);
        }
	}
}();

/* Ported to JavaScript from Java 07/01/14.
 *
 * Ported from C to Java by Dmitry Skiba [sahn0], 23/02/08.
 * Original: http://cds.xs4all.nl:8081/ecdh/
 */
/* Generic 64-bit integer implementation of Curve25519 ECDH
 * Written by Matthijs van Duin, 200608242056
 * Public domain.
 *
 * Based on work by Daniel J Bernstein, http://cr.yp.to/ecdh.html
 */

var curve25519 = function () {

    //region Constants

    var KEY_SIZE = 32;

    /* array length */
    var UNPACKED_SIZE = 16;

    /* group order (a prime near 2^252+2^124) */
    var ORDER = [
        237, 211, 245, 92,
        26, 99, 18, 88,
        214, 156, 247, 162,
        222, 249, 222, 20,
        0, 0, 0, 0,
        0, 0, 0, 16
    ];

    /* smallest multiple of the order that's >= 2^255 */
    var ORDER_TIMES_8 = [
        104, 159, 174, 231,
        210, 24, 147, 192,
        178, 230, 188, 23,
        245, 206, 247, 166,
        0, 0, 0, 0,
        0, 0, 0, 128
    ];

    /* constants 2Gy and 1/(2Gy) */
    var BASE_2Y = [
        22587, 610, 29883, 44076,
        15515, 9479, 25859, 56197,
        23910, 4462, 17831, 16322,
        62102, 36542, 52412, 16035
    ];

    var BASE_R2Y = [
        5744, 16384, 61977, 54121,
        8776, 18501, 26522, 34893,
        23833, 5823, 55924, 58749,
        24147, 14085, 13606, 6080
    ];

    var C1 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    var C9 = [9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    var C486671 = [0x6D0F, 0x0007, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    var C39420360 = [0x81C8, 0x0259, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

    var P25 = 33554431; /* (1 << 25) - 1 */
    var P26 = 67108863; /* (1 << 26) - 1 */

    //#endregion

    //region Key Agreement

    /* Private key clamping
     *   k [out] your private key for key agreement
     *   k  [in]  32 random bytes
     */
    function clamp (k) {
        k[31] &= 0x7F;
        k[31] |= 0x40;
        k[ 0] &= 0xF8;
    }

    //endregion

    //region radix 2^8 math

    function cpy32 (d, s) {
        for (var i = 0; i < 32; i++)
            d[i] = s[i];
    }

    /* p[m..n+m-1] = q[m..n+m-1] + z * x */
    /* n is the size of x */
    /* n+m is the size of p and q */
    function mula_small (p, q, m, x, n, z) {
        m = m | 0;
        n = n | 0;
        z = z | 0;

        var v = 0;
        for (var i = 0; i < n; ++i) {
            v += (q[i + m] & 0xFF) + z * (x[i] & 0xFF);
            p[i + m] = (v & 0xFF);
            v >>= 8;
        }

        return v;
    }

    /* p += x * y * z  where z is a small integer
     * x is size 32, y is size t, p is size 32+t
     * y is allowed to overlap with p+32 if you don't care about the upper half  */
    function mula32 (p, x, y, t, z) {
        t = t | 0;
        z = z | 0;

        var n = 31;
        var w = 0;
        var i = 0;
        for (; i < t; i++) {
            var zy = z * (y[i] & 0xFF);
            w += mula_small(p, p, i, x, n, zy) + (p[i+n] & 0xFF) + zy * (x[n] & 0xFF);
            p[i + n] = w & 0xFF;
            w >>= 8;
        }
        p[i + n] = (w + (p[i + n] & 0xFF)) & 0xFF;
        return w >> 8;
    }

    /* divide r (size n) by d (size t), returning quotient q and remainder r
     * quotient is size n-t+1, remainder is size t
     * requires t > 0 && d[t-1] !== 0
     * requires that r[-1] and d[-1] are valid memory locations
     * q may overlap with r+t */
    function divmod (q, r, n, d, t) {
        n = n | 0;
        t = t | 0;

        var rn = 0;
        var dt = (d[t - 1] & 0xFF) << 8;
        if (t > 1)
            dt |= (d[t - 2] & 0xFF);

        while (n-- >= t) {
            var z = (rn << 16) | ((r[n] & 0xFF) << 8);
            if (n > 0)
                z |= (r[n - 1] & 0xFF);

            var i = n - t + 1;
            z /= dt;
            rn += mula_small(r, r, i, d, t, -z);
            q[i] = (z + rn) & 0xFF;
            /* rn is 0 or -1 (underflow) */
            mula_small(r, r, i, d, t, -rn);
            rn = r[n] & 0xFF;
            r[n] = 0;
        }

        r[t-1] = rn & 0xFF;
    }

    function numsize (x, n) {
        while (n-- !== 0 && x[n] === 0) { }
        return n + 1;
    }

    /* Returns x if a contains the gcd, y if b.
     * Also, the returned buffer contains the inverse of a mod b,
     * as 32-byte signed.
     * x and y must have 64 bytes space for temporary use.
     * requires that a[-1] and b[-1] are valid memory locations  */
    function egcd32 (x, y, a, b) {
        var an, bn = 32, qn, i;
        for (i = 0; i < 32; i++)
            x[i] = y[i] = 0;
        x[0] = 1;
        an = numsize(a, 32);
        if (an === 0)
            return y; /* division by zero */
        var temp = new Array(32);
        while (true) {
            qn = bn - an + 1;
            divmod(temp, b, bn, a, an);
            bn = numsize(b, bn);
            if (bn === 0)
                return x;
            mula32(y, x, temp, qn, -1);

            qn = an - bn + 1;
            divmod(temp, a, an, b, bn);
            an = numsize(a, an);
            if (an === 0)
                return y;
            mula32(x, y, temp, qn, -1);
        }
    }

    //endregion

    //region radix 2^25.5 GF(2^255-19) math

    //region pack / unpack

    /* Convert to internal format from little-endian byte format */
    function unpack (x, m) {
        for (var i = 0; i < KEY_SIZE; i += 2)
            x[i / 2] = m[i] & 0xFF | ((m[i + 1] & 0xFF) << 8);
    }

    /* Check if reduced-form input >= 2^255-19 */
    function is_overflow (x) {
        return (
            ((x[0] > P26 - 19)) &&
                ((x[1] & x[3] & x[5] & x[7] & x[9]) === P25) &&
                ((x[2] & x[4] & x[6] & x[8]) === P26)
            ) || (x[9] > P25);
    }

    /* Convert from internal format to little-endian byte format.  The
     * number must be in a reduced form which is output by the following ops:
     *     unpack, mul, sqr
     *     set --  if input in range 0 .. P25
     * If you're unsure if the number is reduced, first multiply it by 1.  */
    function pack (x, m) {
        for (var i = 0; i < UNPACKED_SIZE; ++i) {
            m[2 * i] = x[i] & 0x00FF;
            m[2 * i + 1] = (x[i] & 0xFF00) >> 8;
        }
    }

    //endregion

    function createUnpackedArray () {
        return new Uint16Array(UNPACKED_SIZE);
    }

    /* Copy a number */
    function cpy (d, s) {
        for (var i = 0; i < UNPACKED_SIZE; ++i)
            d[i] = s[i];
    }

    /* Set a number to value, which must be in range -185861411 .. 185861411 */
    function set (d, s) {
        d[0] = s;
        for (var i = 1; i < UNPACKED_SIZE; ++i)
            d[i] = 0;
    }

    /* Add/subtract two numbers.  The inputs must be in reduced form, and the
     * output isn't, so to do another addition or subtraction on the output,
     * first multiply it by one to reduce it. */
    var add = c255laddmodp;
    var sub = c255lsubmodp;

    /* Multiply a number by a small integer in range -185861411 .. 185861411.
     * The output is in reduced form, the input x need not be.  x and xy may point
     * to the same buffer. */
    var mul_small = c255lmulasmall;

    /* Multiply two numbers.  The output is in reduced form, the inputs need not be. */
    var mul = c255lmulmodp;

    /* Square a number.  Optimization of  mul25519(x2, x, x)  */
    var sqr = c255lsqrmodp;

    /* Calculates a reciprocal.  The output is in reduced form, the inputs need not
     * be.  Simply calculates  y = x^(p-2)  so it's not too fast. */
    /* When sqrtassist is true, it instead calculates y = x^((p-5)/8) */
    function recip (y, x, sqrtassist) {
        var t0 = createUnpackedArray();
        var t1 = createUnpackedArray();
        var t2 = createUnpackedArray();
        var t3 = createUnpackedArray();
        var t4 = createUnpackedArray();

        /* the chain for x^(2^255-21) is straight from djb's implementation */
        var i;
        sqr(t1, x); /*  2 === 2 * 1	*/
        sqr(t2, t1); /*  4 === 2 * 2	*/
        sqr(t0, t2); /*  8 === 2 * 4	*/
        mul(t2, t0, x); /*  9 === 8 + 1	*/
        mul(t0, t2, t1); /* 11 === 9 + 2	*/
        sqr(t1, t0); /* 22 === 2 * 11	*/
        mul(t3, t1, t2); /* 31 === 22 + 9 === 2^5   - 2^0	*/
        sqr(t1, t3); /* 2^6   - 2^1	*/
        sqr(t2, t1); /* 2^7   - 2^2	*/
        sqr(t1, t2); /* 2^8   - 2^3	*/
        sqr(t2, t1); /* 2^9   - 2^4	*/
        sqr(t1, t2); /* 2^10  - 2^5	*/
        mul(t2, t1, t3); /* 2^10  - 2^0	*/
        sqr(t1, t2); /* 2^11  - 2^1	*/
        sqr(t3, t1); /* 2^12  - 2^2	*/
        for (i = 1; i < 5; i++) {
            sqr(t1, t3);
            sqr(t3, t1);
        } /* t3 */ /* 2^20  - 2^10	*/
        mul(t1, t3, t2); /* 2^20  - 2^0	*/
        sqr(t3, t1); /* 2^21  - 2^1	*/
        sqr(t4, t3); /* 2^22  - 2^2	*/
        for (i = 1; i < 10; i++) {
            sqr(t3, t4);
            sqr(t4, t3);
        } /* t4 */ /* 2^40  - 2^20	*/
        mul(t3, t4, t1); /* 2^40  - 2^0	*/
        for (i = 0; i < 5; i++) {
            sqr(t1, t3);
            sqr(t3, t1);
        } /* t3 */ /* 2^50  - 2^10	*/
        mul(t1, t3, t2); /* 2^50  - 2^0	*/
        sqr(t2, t1); /* 2^51  - 2^1	*/
        sqr(t3, t2); /* 2^52  - 2^2	*/
        for (i = 1; i < 25; i++) {
            sqr(t2, t3);
            sqr(t3, t2);
        } /* t3 */ /* 2^100 - 2^50 */
        mul(t2, t3, t1); /* 2^100 - 2^0	*/
        sqr(t3, t2); /* 2^101 - 2^1	*/
        sqr(t4, t3); /* 2^102 - 2^2	*/
        for (i = 1; i < 50; i++) {
            sqr(t3, t4);
            sqr(t4, t3);
        } /* t4 */ /* 2^200 - 2^100 */
        mul(t3, t4, t2); /* 2^200 - 2^0	*/
        for (i = 0; i < 25; i++) {
            sqr(t4, t3);
            sqr(t3, t4);
        } /* t3 */ /* 2^250 - 2^50	*/
        mul(t2, t3, t1); /* 2^250 - 2^0	*/
        sqr(t1, t2); /* 2^251 - 2^1	*/
        sqr(t2, t1); /* 2^252 - 2^2	*/
        if (sqrtassist !== 0) {
            mul(y, x, t2); /* 2^252 - 3 */
        } else {
            sqr(t1, t2); /* 2^253 - 2^3	*/
            sqr(t2, t1); /* 2^254 - 2^4	*/
            sqr(t1, t2); /* 2^255 - 2^5	*/
            mul(y, t1, t0); /* 2^255 - 21	*/
        }
    }

    /* checks if x is "negative", requires reduced input */
    function is_negative (x) {
        var isOverflowOrNegative = is_overflow(x) || x[9] < 0;
        var leastSignificantBit = x[0] & 1;
        return ((isOverflowOrNegative ? 1 : 0) ^ leastSignificantBit) & 0xFFFFFFFF;
    }

    /* a square root */
    function sqrt (x, u) {
        var v = createUnpackedArray();
        var t1 = createUnpackedArray();
        var t2 = createUnpackedArray();

        add(t1, u, u); /* t1 = 2u		*/
        recip(v, t1, 1); /* v = (2u)^((p-5)/8)	*/
        sqr(x, v); /* x = v^2		*/
        mul(t2, t1, x); /* t2 = 2uv^2		*/
        sub(t2, t2, C1); /* t2 = 2uv^2-1		*/
        mul(t1, v, t2); /* t1 = v(2uv^2-1)	*/
        mul(x, u, t1); /* x = uv(2uv^2-1)	*/
    }

    //endregion

    //region JavaScript Fast Math

    function c255lsqr8h (a7, a6, a5, a4, a3, a2, a1, a0) {
        var r = [];
        var v;
        r[0] = (v = a0*a0) & 0xFFFF;
        r[1] = (v = ((v / 0x10000) | 0) + 2*a0*a1) & 0xFFFF;
        r[2] = (v = ((v / 0x10000) | 0) + 2*a0*a2 + a1*a1) & 0xFFFF;
        r[3] = (v = ((v / 0x10000) | 0) + 2*a0*a3 + 2*a1*a2) & 0xFFFF;
        r[4] = (v = ((v / 0x10000) | 0) + 2*a0*a4 + 2*a1*a3 + a2*a2) & 0xFFFF;
        r[5] = (v = ((v / 0x10000) | 0) + 2*a0*a5 + 2*a1*a4 + 2*a2*a3) & 0xFFFF;
        r[6] = (v = ((v / 0x10000) | 0) + 2*a0*a6 + 2*a1*a5 + 2*a2*a4 + a3*a3) & 0xFFFF;
        r[7] = (v = ((v / 0x10000) | 0) + 2*a0*a7 + 2*a1*a6 + 2*a2*a5 + 2*a3*a4) & 0xFFFF;
        r[8] = (v = ((v / 0x10000) | 0) + 2*a1*a7 + 2*a2*a6 + 2*a3*a5 + a4*a4) & 0xFFFF;
        r[9] = (v = ((v / 0x10000) | 0) + 2*a2*a7 + 2*a3*a6 + 2*a4*a5) & 0xFFFF;
        r[10] = (v = ((v / 0x10000) | 0) + 2*a3*a7 + 2*a4*a6 + a5*a5) & 0xFFFF;
        r[11] = (v = ((v / 0x10000) | 0) + 2*a4*a7 + 2*a5*a6) & 0xFFFF;
        r[12] = (v = ((v / 0x10000) | 0) + 2*a5*a7 + a6*a6) & 0xFFFF;
        r[13] = (v = ((v / 0x10000) | 0) + 2*a6*a7) & 0xFFFF;
        r[14] = (v = ((v / 0x10000) | 0) + a7*a7) & 0xFFFF;
        r[15] = ((v / 0x10000) | 0);
        return r;
    }

    function c255lsqrmodp (r, a) {
        var x = c255lsqr8h(a[15], a[14], a[13], a[12], a[11], a[10], a[9], a[8]);
        var z = c255lsqr8h(a[7], a[6], a[5], a[4], a[3], a[2], a[1], a[0]);
        var y = c255lsqr8h(a[15] + a[7], a[14] + a[6], a[13] + a[5], a[12] + a[4], a[11] + a[3], a[10] + a[2], a[9] + a[1], a[8] + a[0]);

        var v;
        r[0] = (v = 0x800000 + z[0] + (y[8] -x[8] -z[8] + x[0] -0x80) * 38) & 0xFFFF;
        r[1] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[1] + (y[9] -x[9] -z[9] + x[1]) * 38) & 0xFFFF;
        r[2] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[2] + (y[10] -x[10] -z[10] + x[2]) * 38) & 0xFFFF;
        r[3] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[3] + (y[11] -x[11] -z[11] + x[3]) * 38) & 0xFFFF;
        r[4] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[4] + (y[12] -x[12] -z[12] + x[4]) * 38) & 0xFFFF;
        r[5] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[5] + (y[13] -x[13] -z[13] + x[5]) * 38) & 0xFFFF;
        r[6] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[6] + (y[14] -x[14] -z[14] + x[6]) * 38) & 0xFFFF;
        r[7] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[7] + (y[15] -x[15] -z[15] + x[7]) * 38) & 0xFFFF;
        r[8] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[8] + y[0] -x[0] -z[0] + x[8] * 38) & 0xFFFF;
        r[9] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[9] + y[1] -x[1] -z[1] + x[9] * 38) & 0xFFFF;
        r[10] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[10] + y[2] -x[2] -z[2] + x[10] * 38) & 0xFFFF;
        r[11] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[11] + y[3] -x[3] -z[3] + x[11] * 38) & 0xFFFF;
        r[12] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[12] + y[4] -x[4] -z[4] + x[12] * 38) & 0xFFFF;
        r[13] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[13] + y[5] -x[5] -z[5] + x[13] * 38) & 0xFFFF;
        r[14] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[14] + y[6] -x[6] -z[6] + x[14] * 38) & 0xFFFF;
        var r15 = 0x7fff80 + ((v / 0x10000) | 0) + z[15] + y[7] -x[7] -z[7] + x[15] * 38;
        c255lreduce(r, r15);
    }

    function c255lmul8h (a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0) {
        var r = [];
        var v;
        r[0] = (v = a0*b0) & 0xFFFF;
        r[1] = (v = ((v / 0x10000) | 0) + a0*b1 + a1*b0) & 0xFFFF;
        r[2] = (v = ((v / 0x10000) | 0) + a0*b2 + a1*b1 + a2*b0) & 0xFFFF;
        r[3] = (v = ((v / 0x10000) | 0) + a0*b3 + a1*b2 + a2*b1 + a3*b0) & 0xFFFF;
        r[4] = (v = ((v / 0x10000) | 0) + a0*b4 + a1*b3 + a2*b2 + a3*b1 + a4*b0) & 0xFFFF;
        r[5] = (v = ((v / 0x10000) | 0) + a0*b5 + a1*b4 + a2*b3 + a3*b2 + a4*b1 + a5*b0) & 0xFFFF;
        r[6] = (v = ((v / 0x10000) | 0) + a0*b6 + a1*b5 + a2*b4 + a3*b3 + a4*b2 + a5*b1 + a6*b0) & 0xFFFF;
        r[7] = (v = ((v / 0x10000) | 0) + a0*b7 + a1*b6 + a2*b5 + a3*b4 + a4*b3 + a5*b2 + a6*b1 + a7*b0) & 0xFFFF;
        r[8] = (v = ((v / 0x10000) | 0) + a1*b7 + a2*b6 + a3*b5 + a4*b4 + a5*b3 + a6*b2 + a7*b1) & 0xFFFF;
        r[9] = (v = ((v / 0x10000) | 0) + a2*b7 + a3*b6 + a4*b5 + a5*b4 + a6*b3 + a7*b2) & 0xFFFF;
        r[10] = (v = ((v / 0x10000) | 0) + a3*b7 + a4*b6 + a5*b5 + a6*b4 + a7*b3) & 0xFFFF;
        r[11] = (v = ((v / 0x10000) | 0) + a4*b7 + a5*b6 + a6*b5 + a7*b4) & 0xFFFF;
        r[12] = (v = ((v / 0x10000) | 0) + a5*b7 + a6*b6 + a7*b5) & 0xFFFF;
        r[13] = (v = ((v / 0x10000) | 0) + a6*b7 + a7*b6) & 0xFFFF;
        r[14] = (v = ((v / 0x10000) | 0) + a7*b7) & 0xFFFF;
        r[15] = ((v / 0x10000) | 0);
        return r;
    }

    function c255lmulmodp (r, a, b) {
        // Karatsuba multiplication scheme: x*y = (b^2+b)*x1*y1 - b*(x1-x0)*(y1-y0) + (b+1)*x0*y0
        var x = c255lmul8h(a[15], a[14], a[13], a[12], a[11], a[10], a[9], a[8], b[15], b[14], b[13], b[12], b[11], b[10], b[9], b[8]);
        var z = c255lmul8h(a[7], a[6], a[5], a[4], a[3], a[2], a[1], a[0], b[7], b[6], b[5], b[4], b[3], b[2], b[1], b[0]);
        var y = c255lmul8h(a[15] + a[7], a[14] + a[6], a[13] + a[5], a[12] + a[4], a[11] + a[3], a[10] + a[2], a[9] + a[1], a[8] + a[0],
            b[15] + b[7], b[14] + b[6], b[13] + b[5], b[12] + b[4], b[11] + b[3], b[10] + b[2], b[9] + b[1], b[8] + b[0]);

        var v;
        r[0] = (v = 0x800000 + z[0] + (y[8] -x[8] -z[8] + x[0] -0x80) * 38) & 0xFFFF;
        r[1] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[1] + (y[9] -x[9] -z[9] + x[1]) * 38) & 0xFFFF;
        r[2] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[2] + (y[10] -x[10] -z[10] + x[2]) * 38) & 0xFFFF;
        r[3] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[3] + (y[11] -x[11] -z[11] + x[3]) * 38) & 0xFFFF;
        r[4] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[4] + (y[12] -x[12] -z[12] + x[4]) * 38) & 0xFFFF;
        r[5] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[5] + (y[13] -x[13] -z[13] + x[5]) * 38) & 0xFFFF;
        r[6] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[6] + (y[14] -x[14] -z[14] + x[6]) * 38) & 0xFFFF;
        r[7] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[7] + (y[15] -x[15] -z[15] + x[7]) * 38) & 0xFFFF;
        r[8] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[8] + y[0] -x[0] -z[0] + x[8] * 38) & 0xFFFF;
        r[9] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[9] + y[1] -x[1] -z[1] + x[9] * 38) & 0xFFFF;
        r[10] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[10] + y[2] -x[2] -z[2] + x[10] * 38) & 0xFFFF;
        r[11] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[11] + y[3] -x[3] -z[3] + x[11] * 38) & 0xFFFF;
        r[12] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[12] + y[4] -x[4] -z[4] + x[12] * 38) & 0xFFFF;
        r[13] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[13] + y[5] -x[5] -z[5] + x[13] * 38) & 0xFFFF;
        r[14] = (v = 0x7fff80 + ((v / 0x10000) | 0) + z[14] + y[6] -x[6] -z[6] + x[14] * 38) & 0xFFFF;
        var r15 = 0x7fff80 + ((v / 0x10000) | 0) + z[15] + y[7] -x[7] -z[7] + x[15] * 38;
        c255lreduce(r, r15);
    }

    function c255lreduce (a, a15) {
        var v = a15;
        a[15] = v & 0x7FFF;
        v = ((v / 0x8000) | 0) * 19;
        for (var i = 0; i <= 14; ++i) {
            a[i] = (v += a[i]) & 0xFFFF;
            v = ((v / 0x10000) | 0);
        }

        a[15] += v;
    }

    function c255laddmodp (r, a, b) {
        var v;
        r[0] = (v = (((a[15] / 0x8000) | 0) + ((b[15] / 0x8000) | 0)) * 19 + a[0] + b[0]) & 0xFFFF;
        for (var i = 1; i <= 14; ++i)
            r[i] = (v = ((v / 0x10000) | 0) + a[i] + b[i]) & 0xFFFF;

        r[15] = ((v / 0x10000) | 0) + (a[15] & 0x7FFF) + (b[15] & 0x7FFF);
    }

    function c255lsubmodp (r, a, b) {
        var v;
        r[0] = (v = 0x80000 + (((a[15] / 0x8000) | 0) - ((b[15] / 0x8000) | 0) - 1) * 19 + a[0] - b[0]) & 0xFFFF;
        for (var i = 1; i <= 14; ++i)
            r[i] = (v = ((v / 0x10000) | 0) + 0x7fff8 + a[i] - b[i]) & 0xFFFF;

        r[15] = ((v / 0x10000) | 0) + 0x7ff8 + (a[15] & 0x7FFF) - (b[15] & 0x7FFF);
    }

    function c255lmulasmall (r, a, m) {
        var v;
        r[0] = (v = a[0] * m) & 0xFFFF;
        for (var i = 1; i <= 14; ++i)
            r[i] = (v = ((v / 0x10000) | 0) + a[i]*m) & 0xFFFF;

        var r15 = ((v / 0x10000) | 0) + a[15]*m;
        c255lreduce(r, r15);
    }

    //endregion

    /********************* Elliptic curve *********************/

    /* y^2 = x^3 + 486662 x^2 + x  over GF(2^255-19) */

    /* t1 = ax + az
     * t2 = ax - az  */
    function mont_prep (t1, t2, ax, az) {
        add(t1, ax, az);
        sub(t2, ax, az);
    }

    /* A = P + Q   where
     *  X(A) = ax/az
     *  X(P) = (t1+t2)/(t1-t2)
     *  X(Q) = (t3+t4)/(t3-t4)
     *  X(P-Q) = dx
     * clobbers t1 and t2, preserves t3 and t4  */
    function mont_add (t1, t2, t3, t4, ax, az, dx) {
        mul(ax, t2, t3);
        mul(az, t1, t4);
        add(t1, ax, az);
        sub(t2, ax, az);
        sqr(ax, t1);
        sqr(t1, t2);
        mul(az, t1, dx);
    }

    /* B = 2 * Q   where
     *  X(B) = bx/bz
     *  X(Q) = (t3+t4)/(t3-t4)
     * clobbers t1 and t2, preserves t3 and t4  */
    function mont_dbl (t1, t2, t3, t4, bx, bz) {
        sqr(t1, t3);
        sqr(t2, t4);
        mul(bx, t1, t2);
        sub(t2, t1, t2);
        mul_small(bz, t2, 121665);
        add(t1, t1, bz);
        mul(bz, t1, t2);
    }

    /* Y^2 = X^3 + 486662 X^2 + X
     * t is a temporary  */
    function x_to_y2 (t, y2, x) {
        sqr(t, x);
        mul_small(y2, x, 486662);
        add(t, t, y2);
        add(t, t, C1);
        mul(y2, t, x);
    }

    /* P = kG   and  s = sign(P)/k  */
    function core (Px, s, k, Gx) {
        var dx = createUnpackedArray();
        var t1 = createUnpackedArray();
        var t2 = createUnpackedArray();
        var t3 = createUnpackedArray();
        var t4 = createUnpackedArray();
        var x = [createUnpackedArray(), createUnpackedArray()];
        var z = [createUnpackedArray(), createUnpackedArray()];
        var i, j;

        /* unpack the base */
        if (Gx !== null)
            unpack(dx, Gx);
        else
            set(dx, 9);

        /* 0G = point-at-infinity */
        set(x[0], 1);
        set(z[0], 0);

        /* 1G = G */
        cpy(x[1], dx);
        set(z[1], 1);

        for (i = 32; i-- !== 0;) {
            for (j = 8; j-- !== 0;) {
                /* swap arguments depending on bit */
                var bit1 = (k[i] & 0xFF) >> j & 1;
                var bit0 = ~(k[i] & 0xFF) >> j & 1;
                var ax = x[bit0];
                var az = z[bit0];
                var bx = x[bit1];
                var bz = z[bit1];

                /* a' = a + b	*/
                /* b' = 2 b	*/
                mont_prep(t1, t2, ax, az);
                mont_prep(t3, t4, bx, bz);
                mont_add(t1, t2, t3, t4, ax, az, dx);
                mont_dbl(t1, t2, t3, t4, bx, bz);
            }
        }

        recip(t1, z[0], 0);
        mul(dx, x[0], t1);

        pack(dx, Px);

        /* calculate s such that s abs(P) = G  .. assumes G is std base point */
        if (s !== null) {
            x_to_y2(t2, t1, dx); /* t1 = Py^2  */
            recip(t3, z[1], 0); /* where Q=P+G ... */
            mul(t2, x[1], t3); /* t2 = Qx  */
            add(t2, t2, dx); /* t2 = Qx + Px  */
            add(t2, t2, C486671); /* t2 = Qx + Px + Gx + 486662  */
            sub(dx, dx, C9); /* dx = Px - Gx  */
            sqr(t3, dx); /* t3 = (Px - Gx)^2  */
            mul(dx, t2, t3); /* dx = t2 (Px - Gx)^2  */
            sub(dx, dx, t1); /* dx = t2 (Px - Gx)^2 - Py^2  */
            sub(dx, dx, C39420360); /* dx = t2 (Px - Gx)^2 - Py^2 - Gy^2  */
            mul(t1, dx, BASE_R2Y); /* t1 = -Py  */

            if (is_negative(t1) !== 0)    /* sign is 1, so just copy  */
                cpy32(s, k);
            else            /* sign is -1, so negate  */
                mula_small(s, ORDER_TIMES_8, 0, k, 32, -1);

            /* reduce s mod q
             * (is this needed?  do it just in case, it's fast anyway) */
            //divmod((dstptr) t1, s, 32, order25519, 32);

            /* take reciprocal of s mod q */
            var temp1 = new Array(32);
            var temp2 = new Array(64);
            var temp3 = new Array(64);
            cpy32(temp1, ORDER);
            cpy32(s, egcd32(temp2, temp3, s, temp1));
            if ((s[31] & 0x80) !== 0)
                mula_small(s, s, 0, ORDER, 32, 1);

        }
    }

    /********* DIGITAL SIGNATURES *********/

    /* deterministic EC-KCDSA
     *
     *    s is the private key for signing
     *    P is the corresponding public key
     *    Z is the context data (signer public key or certificate, etc)
     *
     * signing:
     *
     *    m = hash(Z, message)
     *    x = hash(m, s)
     *    keygen25519(Y, NULL, x);
     *    r = hash(Y);
     *    h = m XOR r
     *    sign25519(v, h, x, s);
     *
     *    output (v,r) as the signature
     *
     * verification:
     *
     *    m = hash(Z, message);
     *    h = m XOR r
     *    verify25519(Y, v, h, P)
     *
     *    confirm  r === hash(Y)
     *
     * It would seem to me that it would be simpler to have the signer directly do
     * h = hash(m, Y) and send that to the recipient instead of r, who can verify
     * the signature by checking h === hash(m, Y).  If there are any problems with
     * such a scheme, please let me know.
     *
     * Also, EC-KCDSA (like most DS algorithms) picks x random, which is a waste of
     * perfectly good entropy, but does allow Y to be calculated in advance of (or
     * parallel to) hashing the message.
     */

    /* Signature generation primitive, calculates (x-h)s mod q
     *   h  [in]  signature hash (of message, signature pub key, and context data)
     *   x  [in]  signature private key
     *   s  [in]  private key for signing
     * returns signature value on success, undefined on failure (use different x or h)
     */

    function sign (h, x, s) {
        // v = (x - h) s  mod q
        var w, i;
        var h1 = new Array(32)
        var x1 = new Array(32);
        var tmp1 = new Array(64);
        var tmp2 = new Array(64);

        // Don't clobber the arguments, be nice!
        cpy32(h1, h);
        cpy32(x1, x);

        // Reduce modulo group order
        var tmp3 = new Array(32);
        divmod(tmp3, h1, 32, ORDER, 32);
        divmod(tmp3, x1, 32, ORDER, 32);

        // v = x1 - h1
        // If v is negative, add the group order to it to become positive.
        // If v was already positive we don't have to worry about overflow
        // when adding the order because v < ORDER and 2*ORDER < 2^256
        var v = new Array(32);
        mula_small(v, x1, 0, h1, 32, -1);
        mula_small(v, v , 0, ORDER, 32, 1);

        // tmp1 = (x-h)*s mod q
        mula32(tmp1, v, s, 32, 1);
        divmod(tmp2, tmp1, 64, ORDER, 32);

        for (w = 0, i = 0; i < 32; i++)
            w |= v[i] = tmp1[i];

        return w !== 0 ? v : undefined;
    }

    /* Signature verification primitive, calculates Y = vP + hG
     *   v  [in]  signature value
     *   h  [in]  signature hash
     *   P  [in]  public key
     *   Returns signature public key
     */
    function verify (v, h, P) {
        /* Y = v abs(P) + h G  */
        var d = new Array(32);
        var p = [createUnpackedArray(), createUnpackedArray()];
        var s = [createUnpackedArray(), createUnpackedArray()];
        var yx = [createUnpackedArray(), createUnpackedArray(), createUnpackedArray()];
        var yz = [createUnpackedArray(), createUnpackedArray(), createUnpackedArray()];
        var t1 = [createUnpackedArray(), createUnpackedArray(), createUnpackedArray()];
        var t2 = [createUnpackedArray(), createUnpackedArray(), createUnpackedArray()];

        var vi = 0, hi = 0, di = 0, nvh = 0, i, j, k;

        /* set p[0] to G and p[1] to P  */

        set(p[0], 9);
        unpack(p[1], P);

        /* set s[0] to P+G and s[1] to P-G  */

        /* s[0] = (Py^2 + Gy^2 - 2 Py Gy)/(Px - Gx)^2 - Px - Gx - 486662  */
        /* s[1] = (Py^2 + Gy^2 + 2 Py Gy)/(Px - Gx)^2 - Px - Gx - 486662  */

        x_to_y2(t1[0], t2[0], p[1]); /* t2[0] = Py^2  */
        sqrt(t1[0], t2[0]); /* t1[0] = Py or -Py  */
        j = is_negative(t1[0]); /*      ... check which  */
        add(t2[0], t2[0], C39420360); /* t2[0] = Py^2 + Gy^2  */
        mul(t2[1], BASE_2Y, t1[0]); /* t2[1] = 2 Py Gy or -2 Py Gy  */
        sub(t1[j], t2[0], t2[1]); /* t1[0] = Py^2 + Gy^2 - 2 Py Gy  */
        add(t1[1 - j], t2[0], t2[1]); /* t1[1] = Py^2 + Gy^2 + 2 Py Gy  */
        cpy(t2[0], p[1]); /* t2[0] = Px  */
        sub(t2[0], t2[0], C9); /* t2[0] = Px - Gx  */
        sqr(t2[1], t2[0]); /* t2[1] = (Px - Gx)^2  */
        recip(t2[0], t2[1], 0); /* t2[0] = 1/(Px - Gx)^2  */
        mul(s[0], t1[0], t2[0]); /* s[0] = t1[0]/(Px - Gx)^2  */
        sub(s[0], s[0], p[1]); /* s[0] = t1[0]/(Px - Gx)^2 - Px  */
        sub(s[0], s[0], C486671); /* s[0] = X(P+G)  */
        mul(s[1], t1[1], t2[0]); /* s[1] = t1[1]/(Px - Gx)^2  */
        sub(s[1], s[1], p[1]); /* s[1] = t1[1]/(Px - Gx)^2 - Px  */
        sub(s[1], s[1], C486671); /* s[1] = X(P-G)  */
        mul_small(s[0], s[0], 1); /* reduce s[0] */
        mul_small(s[1], s[1], 1); /* reduce s[1] */

        /* prepare the chain  */
        for (i = 0; i < 32; i++) {
            vi = (vi >> 8) ^ (v[i] & 0xFF) ^ ((v[i] & 0xFF) << 1);
            hi = (hi >> 8) ^ (h[i] & 0xFF) ^ ((h[i] & 0xFF) << 1);
            nvh = ~(vi ^ hi);
            di = (nvh & (di & 0x80) >> 7) ^ vi;
            di ^= nvh & (di & 0x01) << 1;
            di ^= nvh & (di & 0x02) << 1;
            di ^= nvh & (di & 0x04) << 1;
            di ^= nvh & (di & 0x08) << 1;
            di ^= nvh & (di & 0x10) << 1;
            di ^= nvh & (di & 0x20) << 1;
            di ^= nvh & (di & 0x40) << 1;
            d[i] = di & 0xFF;
        }

        di = ((nvh & (di & 0x80) << 1) ^ vi) >> 8;

        /* initialize state */
        set(yx[0], 1);
        cpy(yx[1], p[di]);
        cpy(yx[2], s[0]);
        set(yz[0], 0);
        set(yz[1], 1);
        set(yz[2], 1);

        /* y[0] is (even)P + (even)G
         * y[1] is (even)P + (odd)G  if current d-bit is 0
         * y[1] is (odd)P + (even)G  if current d-bit is 1
         * y[2] is (odd)P + (odd)G
         */

        vi = 0;
        hi = 0;

        /* and go for it! */
        for (i = 32; i-- !== 0;) {
            vi = (vi << 8) | (v[i] & 0xFF);
            hi = (hi << 8) | (h[i] & 0xFF);
            di = (di << 8) | (d[i] & 0xFF);

            for (j = 8; j-- !== 0;) {
                mont_prep(t1[0], t2[0], yx[0], yz[0]);
                mont_prep(t1[1], t2[1], yx[1], yz[1]);
                mont_prep(t1[2], t2[2], yx[2], yz[2]);

                k = ((vi ^ vi >> 1) >> j & 1)
                    + ((hi ^ hi >> 1) >> j & 1);
                mont_dbl(yx[2], yz[2], t1[k], t2[k], yx[0], yz[0]);

                k = (di >> j & 2) ^ ((di >> j & 1) << 1);
                mont_add(t1[1], t2[1], t1[k], t2[k], yx[1], yz[1],
                    p[di >> j & 1]);

                mont_add(t1[2], t2[2], t1[0], t2[0], yx[2], yz[2],
                    s[((vi ^ hi) >> j & 2) >> 1]);
            }
        }

        k = (vi & 1) + (hi & 1);
        recip(t1[0], yz[k], 0);
        mul(t1[1], yx[k], t1[0]);

        var Y = [];
        pack(t1[1], Y);
        return Y;
    }

    /* Key-pair generation
     *   P  [out] your public key
     *   s  [out] your private key for signing
     *   k  [out] your private key for key agreement
     *   k  [in]  32 random bytes
     * s may be NULL if you don't care
     *
     * WARNING: if s is not NULL, this function has data-dependent timing */
    function keygen (k) {
        var P = [];
        var s = [];
        k = k || [];
        clamp(k);
        core(P, s, k, null);

        return { p: P, s: s, k: k };
    }

    return {
        sign: sign,
        verify: verify,
        keygen: keygen
    };
}();

// Copyright (c) 2007 Michele Bini
// Konstantin Welke, 2008:
// - moved into .js file, renamed all c255lname to curve25519_name
// - added curve25519_clamp()
// - functions to read from/to 8bit string
// - removed base32/hex functions (cleanup)
// - removed setbit function (cleanup, had a bug anyway)
// BloodyRookie 2014:
// - ported part of the java implementation by Dmitry Skiba to js and merged into this file 
// - profiled for higher speed
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */
//
// The original curve25519 library was released into the public domain
// by Daniel J. Bernstein

curve25519_zero = function() {
  return [0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0];
}

curve25519_one = function() {
  return [1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0];
}

curve25519_two = function() {
  return [2,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0];
}

curve25519_nine = function() {
  return [9,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0];
}

curve25519_486671 = function() {
  return [27919,7,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0];
}

curve25519_39420360 = function() {
  return [33224,601,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0];
}

curve25519_r2y = function() {
  return [0x1670,0x4000,0xf219,0xd369,0x2248,0x4845,0x679a,0x884d,0x5d19,0x16bf,0xda74,0xe57d,0x5e53,0x3705,0x3526,0x17c0];
}

curve25519_2y = function() {
  return [0x583b,0x0262,0x74bb,0xac2c,0x3c9b,0x2507,0x6503,0xdb85,0x5d66,0x116e,0x45a7,0x3fc2,0xf296,0x8ebe,0xccbc,0x3ea3];
}

curve25519_clamp = function(curve) {
  curve[0] &= 0xFFF8;
  curve[15] &= 0x7FFF;
  curve[15] |= 0x4000;
  return curve;
}

curve25519_getbit = function(curve, c) {
  return ~~(curve[~~(c / 16)] / Math.pow(2, c % 16)) % 2;
}
  
curve25519_prime = [0xffff-18, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0x7fff];

/* group order (a prime near 2^252+2^124) */
curve25519_order = [
	    237, 211, 245, 92, 26, 99, 18, 88, 214, 156, 247, 162, 222, 249, 222, 20,
	    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16];

curve25519_order_times_8 = [
	    104, 159, 174, 231, 210, 24, 147, 192, 178, 230, 188, 23, 245, 206, 247, 166,
	    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128];

curve25519_convertToByteArray = function(a) {
	var b = new Int8Array(32);
	var i;
	for (i=0; i<16; i++)
	{
		b[2*i] = a[i] & 0xff;
		b[2*i+1] = a[i] >> 8;
	}
	
	return b;
}

curve25519_convertToShortArray = function(a) {
	var b = new Array(16);
	var i, val1, val2;
	for (i=0; i<16; i++)
	{
		val1 = a[i*2];
		if (val1 < 0)
		{
			val1 +=256;
		}
		val2 = a[i*2+1];
		if (val2 < 0)
		{
			val2 +=256;
		}
		b[i] = val1 + val2*256;
	}
	return b;
	
}

curve25519_fillShortArray = function(src, dest) {
	var i;
	for (i=0; i<16; i++)
	{
		dest[i] = src[i];
	}
}

curve25519_fillByteArray = function(src, dest) {
	var i;
	for (i=0; i<32; i++)
	{
		dest[i] = src[i];
	}
}

curve25519_log16 = function(text, a) {
	var b = shortArray_to_hex_string(a);
	addText(text + b);
}

curve25519_log32 = function(text, a) {
	var b = byteArray_to_hex_string(a);
	addText(text + b);
}

curve25519_cpy32 = function(a) {
	var b = new Int8Array(32);
	for (i = 0; i < 32; i++)
	{
		b[i] = a[i];
	}
	return b;
}

curve25519_mula_small = function(p, q, m, x, n, z) {
	var v=0;
	for (j=0; j<n; ++j) 
	{
		v += (q[j+m] & 0xFF) + z * (x[j] & 0xFF);
		p[j+m] = (v & 0xFF);
		v >>= 8;
	}
	return v;		
}

curve25519_mula32 = function(p, x, y, t, z) {
	var n = 31;
	var w = 0;
	for (i=0; i < t; i++) 
	{
		zy = z * (y[i] & 0xFF);
		w += curve25519_mula_small(p, p, i, x, n, zy) + (p[i+n] & 0xFF) + zy * (x[n] & 0xFF);
		p[i+n] = (w & 0xFF);
		w >>= 8;
	}
	p[i+n] = ((w + (p[i+n] & 0xFF)) & 0xFF);
	return w >> 8;
}

curve25519_divmod = function(q, r, n, d, t) {
	var rn = 0, z=0;
	var dt = ((d[t-1] & 0xFF) << 8);
	if (t>1) 
	{
		dt |= (d[t-2] & 0xFF);
	}
	while (n-- >= t) 
	{
		z = (rn << 16) | ((r[n] & 0xFF) << 8);
		if (n>0) 
		{
			z |= (r[n-1] & 0xFF);
		}
		z = parseInt(z/dt);
		rn += curve25519_mula_small(r,r, n-t+1, d, t, -z);
		q[n-t+1] = ((z + rn) & 0xFF); // rn is 0 or -1 (underflow)
		curve25519_mula_small(r,r, n-t+1, d, t, -rn);
		rn = (r[n] & 0xFF);
		r[n] = 0;
	}
	r[t-1] = (rn & 0xFF);
}

curve25519_numsize = function(x, n)  {
	while (n--!=0 && x[n]==0)
		;
	return n+1;
}

curve25519_egcd32 = function(x, y, a, b) {
	var an = 0, bn = 32, qn=0, i=0;
	for (i = 0; i < 32; i++)
	{
		x[i] = y[i] = 0;
	}
	x[0] = 1;
	an = curve25519_numsize(a, 32);
	if (an==0)
	{
		return y;	// division by zero
	}
	temp=new Int8Array(32);
	while (true) 
	{
		qn = bn - an + 1;
		curve25519_divmod(temp, b, bn, a, an);
		bn = curve25519_numsize(b, bn);
		if (bn==0)
		{
			return x;
		}
		curve25519_mula32(y, x, temp, qn, -1);

		qn = an - bn + 1;
		curve25519_divmod(temp, a, an, b, bn);
		an = curve25519_numsize(a, an);
		if (an==0)
		{
			return y;
		}
		curve25519_mula32(x, y, temp, qn, -1);
	}
}

curve25519_compare = function (a ,b) {
  var c;
  for (c = 15; c >= 0; c--) {
    var x = a[c];
    var y = b[c];
    if (x > y) {
      return 1;
    }
    if (x < y) {
      return -1;
    }
  }
  return 0;
}

curve25519_cpy16 = function(a) {
	var r = new Array(16);
	var i;
	for (i=0; i<16;i++)
	{
		r[i] = a[i];
	}
	return r;
}

/***
 * BloodyRookie: odd numbers are negativ
 */
curve25519_isNegative = function(x) {
	return (x[0] & 1);
}

curve25519_isOverflow = function(x) {
	if (x[15] >= 0x8000) return 1;
	if (x[0] >= 0x10000)
	{
		var i;
		for (i=1; i<15; i++)
		{
			if (x[i] < 0xFFFF)
			{
				return 0;
			}
		}
		return 1;
	}
	else	
	{
		return 0;
	}
}

curve25519_sqr8h = function (r, a7, a6, a5, a4, a3, a2, a1, a0) {
  var v=0;
  r[0] = (v = a0*a0) & 0xffff;
  r[1] = (v = ~~(v / 0x10000) + 2*a0*a1) & 0xffff;
  r[2] = (v = ~~(v / 0x10000) + 2*a0*a2 + a1*a1) & 0xffff;
  r[3] = (v = ~~(v / 0x10000) + 2*a0*a3 + 2*a1*a2) & 0xffff;
  r[4] = (v = ~~(v / 0x10000) + 2*a0*a4 + 2*a1*a3 + a2*a2) & 0xffff;
  r[5] = (v = ~~(v / 0x10000) + 2*a0*a5 + 2*a1*a4 + 2*a2*a3) & 0xffff;
  r[6] = (v = ~~(v / 0x10000) + 2*a0*a6 + 2*a1*a5 + 2*a2*a4 + a3*a3) & 0xffff;
  r[7] = (v = ~~(v / 0x10000) + 2*a0*a7 + 2*a1*a6 + 2*a2*a5 + 2*a3*a4) & 0xffff;
  r[8] = (v = ~~(v / 0x10000) + 2*a1*a7 + 2*a2*a6 + 2*a3*a5 + a4*a4) & 0xffff;
  r[9] = (v = ~~(v / 0x10000) + 2*a2*a7 + 2*a3*a6 + 2*a4*a5) & 0xffff;
  r[10] = (v = ~~(v / 0x10000) + 2*a3*a7 + 2*a4*a6 + a5*a5) & 0xffff;
  r[11] = (v = ~~(v / 0x10000) + 2*a4*a7 + 2*a5*a6) & 0xffff;
  r[12] = (v = ~~(v / 0x10000) + 2*a5*a7 + a6*a6) & 0xffff;
  r[13] = (v = ~~(v / 0x10000) + 2*a6*a7) & 0xffff;
  r[14] = (v = ~~(v / 0x10000) + a7*a7) & 0xffff;
  r[15] = ~~(v / 0x10000);
}

curve25519_sqrmodp = function(r, a) {
  var x = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var y = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var z = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  curve25519_sqr8h(x, a[15], a[14], a[13], a[12], a[11], a[10], a[9], a[8]);
  curve25519_sqr8h(z, a[7], a[6], a[5], a[4], a[3], a[2], a[1], a[0]);
  curve25519_sqr8h(y, a[15] + a[7], a[14] + a[6], a[13] + a[5], a[12] + a[4], a[11] + a[3], a[10] + a[2], a[9] + a[1], a[8] + a[0]);
  var v=0;
  r[0] = (v = 0x800000 + z[0] + (y[8] -x[8] -z[8] + x[0] -0x80) * 38) & 0xffff;
  r[1] = (v = 0x7fff80 + ~~(v / 0x10000) + z[1] + (y[9] -x[9] -z[9] + x[1]) * 38) & 0xffff;
  r[2] = (v = 0x7fff80 + ~~(v / 0x10000) + z[2] + (y[10] -x[10] -z[10] + x[2]) * 38) & 0xffff;
  r[3] = (v = 0x7fff80 + ~~(v / 0x10000) + z[3] + (y[11] -x[11] -z[11] + x[3]) * 38) & 0xffff;
  r[4] = (v = 0x7fff80 + ~~(v / 0x10000) + z[4] + (y[12] -x[12] -z[12] + x[4]) * 38) & 0xffff;
  r[5] = (v = 0x7fff80 + ~~(v / 0x10000) + z[5] + (y[13] -x[13] -z[13] + x[5]) * 38) & 0xffff;
  r[6] = (v = 0x7fff80 + ~~(v / 0x10000) + z[6] + (y[14] -x[14] -z[14] + x[6]) * 38) & 0xffff;
  r[7] = (v = 0x7fff80 + ~~(v / 0x10000) + z[7] + (y[15] -x[15] -z[15] + x[7]) * 38) & 0xffff;
  r[8] = (v = 0x7fff80 + ~~(v / 0x10000) + z[8] + y[0] -x[0] -z[0] + x[8] * 38) & 0xffff;
  r[9] = (v = 0x7fff80 + ~~(v / 0x10000) + z[9] + y[1] -x[1] -z[1] + x[9] * 38) & 0xffff;
  r[10] = (v = 0x7fff80 + ~~(v / 0x10000) + z[10] + y[2] -x[2] -z[2] + x[10] * 38) & 0xffff;
  r[11] = (v = 0x7fff80 + ~~(v / 0x10000) + z[11] + y[3] -x[3] -z[3] + x[11] * 38) & 0xffff;
  r[12] = (v = 0x7fff80 + ~~(v / 0x10000) + z[12] + y[4] -x[4] -z[4] + x[12] * 38) & 0xffff;
  r[13] = (v = 0x7fff80 + ~~(v / 0x10000) + z[13] + y[5] -x[5] -z[5] + x[13] * 38) & 0xffff;
  r[14] = (v = 0x7fff80 + ~~(v / 0x10000) + z[14] + y[6] -x[6] -z[6] + x[14] * 38) & 0xffff;
  r[15] = 0x7fff80 + ~~(v / 0x10000) + z[15] + y[7] -x[7] -z[7] + x[15] * 38;
  curve25519_reduce(r);
}

curve25519_mul8h = function(r, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0) {
  var v=0;
  r[0] = (v = a0*b0) & 0xffff;
  r[1] = (v = ~~(v / 0x10000) + a0*b1 + a1*b0) & 0xffff;
  r[2] = (v = ~~(v / 0x10000) + a0*b2 + a1*b1 + a2*b0) & 0xffff;
  r[3] = (v = ~~(v / 0x10000) + a0*b3 + a1*b2 + a2*b1 + a3*b0) & 0xffff;
  r[4] = (v = ~~(v / 0x10000) + a0*b4 + a1*b3 + a2*b2 + a3*b1 + a4*b0) & 0xffff;
  r[5] = (v = ~~(v / 0x10000) + a0*b5 + a1*b4 + a2*b3 + a3*b2 + a4*b1 + a5*b0) & 0xffff;
  r[6] = (v = ~~(v / 0x10000) + a0*b6 + a1*b5 + a2*b4 + a3*b3 + a4*b2 + a5*b1 + a6*b0) & 0xffff;
  r[7] = (v = ~~(v / 0x10000) + a0*b7 + a1*b6 + a2*b5 + a3*b4 + a4*b3 + a5*b2 + a6*b1 + a7*b0) & 0xffff;
  r[8] = (v = ~~(v / 0x10000) + a1*b7 + a2*b6 + a3*b5 + a4*b4 + a5*b3 + a6*b2 + a7*b1) & 0xffff;
  r[9] = (v = ~~(v / 0x10000) + a2*b7 + a3*b6 + a4*b5 + a5*b4 + a6*b3 + a7*b2) & 0xffff;
  r[10] = (v = ~~(v / 0x10000) + a3*b7 + a4*b6 + a5*b5 + a6*b4 + a7*b3) & 0xffff;
  r[11] = (v = ~~(v / 0x10000) + a4*b7 + a5*b6 + a6*b5 + a7*b4) & 0xffff;
  r[12] = (v = ~~(v / 0x10000) + a5*b7 + a6*b6 + a7*b5) & 0xffff;
  r[13] = (v = ~~(v / 0x10000) + a6*b7 + a7*b6) & 0xffff;
  r[14] = (v = ~~(v / 0x10000) + a7*b7) & 0xffff;
  r[15] = ~~(v / 0x10000);
}

curve25519_mulmodp = function(r, a, b) {
  var x = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var y = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var z = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  curve25519_mul8h(x, a[15], a[14], a[13], a[12], a[11], a[10], a[9], a[8], b[15], b[14], b[13], b[12], b[11], b[10], b[9], b[8]);
  curve25519_mul8h(z, a[7], a[6], a[5], a[4], a[3], a[2], a[1], a[0], b[7], b[6], b[5], b[4], b[3], b[2], b[1], b[0]);
  curve25519_mul8h(y, a[15] + a[7], a[14] + a[6], a[13] + a[5], a[12] + a[4], a[11] + a[3], a[10] + a[2], a[9] + a[1], a[8] + a[0],
                   b[15] + b[7], b[14] + b[6], b[13] + b[5], b[12] + b[4], b[11] + b[3], b[10] + b[2], b[9] + b[1], b[8] + b[0]);
  var v=0;
  r[0] = (v = 0x800000 + z[0] + (y[8] -x[8] -z[8] + x[0] -0x80) * 38) & 0xffff;
  r[1] = (v = 0x7fff80 + ~~(v / 0x10000) + z[1] + (y[9] -x[9] -z[9] + x[1]) * 38) & 0xffff;
  r[2] = (v = 0x7fff80 + ~~(v / 0x10000) + z[2] + (y[10] -x[10] -z[10] + x[2]) * 38) & 0xffff;
  r[3] = (v = 0x7fff80 + ~~(v / 0x10000) + z[3] + (y[11] -x[11] -z[11] + x[3]) * 38) & 0xffff;
  r[4] = (v = 0x7fff80 + ~~(v / 0x10000) + z[4] + (y[12] -x[12] -z[12] + x[4]) * 38) & 0xffff;
  r[5] = (v = 0x7fff80 + ~~(v / 0x10000) + z[5] + (y[13] -x[13] -z[13] + x[5]) * 38) & 0xffff;
  r[6] = (v = 0x7fff80 + ~~(v / 0x10000) + z[6] + (y[14] -x[14] -z[14] + x[6]) * 38) & 0xffff;
  r[7] = (v = 0x7fff80 + ~~(v / 0x10000) + z[7] + (y[15] -x[15] -z[15] + x[7]) * 38) & 0xffff;
  r[8] = (v = 0x7fff80 + ~~(v / 0x10000) + z[8] + y[0] -x[0] -z[0] + x[8] * 38) & 0xffff;
  r[9] = (v = 0x7fff80 + ~~(v / 0x10000) + z[9] + y[1] -x[1] -z[1] + x[9] * 38) & 0xffff;
  r[10] = (v = 0x7fff80 + ~~(v / 0x10000) + z[10] + y[2] -x[2] -z[2] + x[10] * 38) & 0xffff;
  r[11] = (v = 0x7fff80 + ~~(v / 0x10000) + z[11] + y[3] -x[3] -z[3] + x[11] * 38) & 0xffff;
  r[12] = (v = 0x7fff80 + ~~(v / 0x10000) + z[12] + y[4] -x[4] -z[4] + x[12] * 38) & 0xffff;
  r[13] = (v = 0x7fff80 + ~~(v / 0x10000) + z[13] + y[5] -x[5] -z[5] + x[13] * 38) & 0xffff;
  r[14] = (v = 0x7fff80 + ~~(v / 0x10000) + z[14] + y[6] -x[6] -z[6] + x[14] * 38) & 0xffff;
  r[15] = 0x7fff80 + ~~(v / 0x10000) + z[15] + y[7] -x[7] -z[7] + x[15] * 38;
  curve25519_reduce(r);
}

curve25519_mulasmall = function(r, a, m) {
  var v=0;
  r[0] = (v = a[0] * m) & 0xffff;
  r[1] = (v = ~~(v / 0x10000) + a[1]*m) & 0xffff;
  r[2] = (v = ~~(v / 0x10000) + a[2]*m) & 0xffff;
  r[3] = (v = ~~(v / 0x10000) + a[3]*m) & 0xffff;
  r[4] = (v = ~~(v / 0x10000) + a[4]*m) & 0xffff;
  r[5] = (v = ~~(v / 0x10000) + a[5]*m) & 0xffff;
  r[6] = (v = ~~(v / 0x10000) + a[6]*m) & 0xffff;
  r[7] = (v = ~~(v / 0x10000) + a[7]*m) & 0xffff;
  r[8] = (v = ~~(v / 0x10000) + a[8]*m) & 0xffff;
  r[9] = (v = ~~(v / 0x10000) + a[9]*m) & 0xffff;
  r[10] = (v = ~~(v / 0x10000) + a[10]*m) & 0xffff;
  r[11] = (v = ~~(v / 0x10000) + a[11]*m) & 0xffff;
  r[12] = (v = ~~(v / 0x10000) + a[12]*m) & 0xffff;
  r[13] = (v = ~~(v / 0x10000) + a[13]*m) & 0xffff;
  r[14] = (v = ~~(v / 0x10000) + a[14]*m) & 0xffff;
  r[15] = ~~(v / 0x10000) + a[15]*m;
  curve25519_reduce(r);
}

curve25519_addmodp = function(r, a, b) {
  var v=0;
  r[0] = (v = (~~(a[15] / 0x8000) + ~~(b[15] / 0x8000)) * 19 + a[0] + b[0]) & 0xffff;
  r[1] = (v = ~~(v / 0x10000) + a[1] + b[1]) & 0xffff;
  r[2] = (v = ~~(v / 0x10000) + a[2] + b[2]) & 0xffff;
  r[3] = (v = ~~(v / 0x10000) + a[3] + b[3]) & 0xffff;
  r[4] = (v = ~~(v / 0x10000) + a[4] + b[4]) & 0xffff;
  r[5] = (v = ~~(v / 0x10000) + a[5] + b[5]) & 0xffff;
  r[6] = (v = ~~(v / 0x10000) + a[6] + b[6]) & 0xffff;
  r[7] = (v = ~~(v / 0x10000) + a[7] + b[7]) & 0xffff;
  r[8] = (v = ~~(v / 0x10000) + a[8] + b[8]) & 0xffff;
  r[9] = (v = ~~(v / 0x10000) + a[9] + b[9]) & 0xffff;
  r[10] = (v = ~~(v / 0x10000) + a[10] + b[10]) & 0xffff;
  r[11] = (v = ~~(v / 0x10000) + a[11] + b[11]) & 0xffff;
  r[12] = (v = ~~(v / 0x10000) + a[12] + b[12]) & 0xffff;
  r[13] = (v = ~~(v / 0x10000) + a[13] + b[13]) & 0xffff;
  r[14] = (v = ~~(v / 0x10000) + a[14] + b[14]) & 0xffff;
  r[15] = ~~(v / 0x10000) + a[15] % 0x8000 + b[15] % 0x8000;
}

curve25519_submodp = function(r, a, b) {
  var v=0;
  r[0] = (v = 0x80000 + (~~(a[15] / 0x8000) - ~~(b[15] / 0x8000) - 1) * 19 + a[0] - b[0]) & 0xffff;
  r[1] = (v = ~~(v / 0x10000) + 0x7fff8 + a[1] - b[1]) & 0xffff;
  r[2] = (v = ~~(v / 0x10000) + 0x7fff8 + a[2] - b[2]) & 0xffff;
  r[3] = (v = ~~(v / 0x10000) + 0x7fff8 + a[3] - b[3]) & 0xffff;
  r[4] = (v = ~~(v / 0x10000) + 0x7fff8 + a[4] - b[4]) & 0xffff;
  r[5] = (v = ~~(v / 0x10000) + 0x7fff8 + a[5] - b[5]) & 0xffff;
  r[6] = (v = ~~(v / 0x10000) + 0x7fff8 + a[6] - b[6]) & 0xffff;
  r[7] = (v = ~~(v / 0x10000) + 0x7fff8 + a[7] - b[7]) & 0xffff;
  r[8] = (v = ~~(v / 0x10000) + 0x7fff8 + a[8] - b[8]) & 0xffff;
  r[9] = (v = ~~(v / 0x10000) + 0x7fff8 + a[9] - b[9]) & 0xffff;
  r[10] = (v = ~~(v / 0x10000) + 0x7fff8 + a[10] - b[10]) & 0xffff;
  r[11] = (v = ~~(v / 0x10000) + 0x7fff8 + a[11] - b[11]) & 0xffff;
  r[12] = (v = ~~(v / 0x10000) + 0x7fff8 + a[12] - b[12]) & 0xffff;
  r[13] = (v = ~~(v / 0x10000) + 0x7fff8 + a[13] - b[13]) & 0xffff;
  r[14] = (v = ~~(v / 0x10000) + 0x7fff8 + a[14] - b[14]) & 0xffff;
  r[15] = ~~(v / 0x10000) + 0x7ff8 + a[15]%0x8000 - b[15]%0x8000;
}
/****
 * BloodyRookie: a^-1 is found via Fermats little theorem:
 * a^p congruent a mod p and therefore a^(p-2) congruent a^-1 mod p
 */
curve25519_invmodp = function (r, a, sqrtassist) {
  var r1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r3 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r4 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r5 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var i=0;
  curve25519_sqrmodp(r2, a);					//  2 == 2 * 1	
  curve25519_sqrmodp(r3, r2);					//  4 == 2 * 2	
  curve25519_sqrmodp(r1, r3);					//  8 == 2 * 4	
  curve25519_mulmodp(r3, r1, a);				//  9 == 8 + 1	
  curve25519_mulmodp(r1, r3, r2);				// 11 == 9 + 2	
  curve25519_sqrmodp(r2, r1);					// 22 == 2 * 11	
  curve25519_mulmodp(r4, r2, r3);				// 31 == 22 + 9
 												//	== 2^5   - 2^0	
  curve25519_sqrmodp(r2, r4);					// 2^6   - 2^1	
  curve25519_sqrmodp(r3, r2);					// 2^7   - 2^2	
  curve25519_sqrmodp(r2, r3);					// 2^8   - 2^3	
  curve25519_sqrmodp(r3, r2);					// 2^9   - 2^4	
  curve25519_sqrmodp(r2, r3);					// 2^10  - 2^5	
  curve25519_mulmodp(r3, r2, r4);				// 2^10  - 2^0	
  curve25519_sqrmodp(r2, r3);					// 2^11  - 2^1	
  curve25519_sqrmodp(r4, r2);					// 2^12  - 2^2	
  for (i = 1; i < 5; i++) {
	curve25519_sqrmodp(r2, r4);
	curve25519_sqrmodp(r4, r2);
  } 											// 2^20  - 2^10	
  curve25519_mulmodp(r2, r4, r3);				// 2^20  - 2^0	
  curve25519_sqrmodp(r4, r2);					// 2^21  - 2^1	
  curve25519_sqrmodp(r5, r4);					// 2^22  - 2^2	
  for (i = 1; i < 10; i++) {
	curve25519_sqrmodp(r4, r5);
	curve25519_sqrmodp(r5, r4);
  } 											// 2^40  - 2^20	
  curve25519_mulmodp(r4, r5, r2);				// 2^40  - 2^0	
  for (i = 0; i < 5; i++) {
	curve25519_sqrmodp(r2, r4);
	curve25519_sqrmodp(r4, r2);
  } 											// 2^50  - 2^10	
  curve25519_mulmodp(r2, r4, r3);				// 2^50  - 2^0	
  curve25519_sqrmodp(r3, r2);					// 2^51  - 2^1	
  curve25519_sqrmodp(r4, r3);					// 2^52  - 2^2	
  for (i = 1; i < 25; i++) {
	curve25519_sqrmodp(r3, r4);
	curve25519_sqrmodp(r4, r3);
  } 											// 2^100 - 2^50 
  curve25519_mulmodp(r3, r4, r2);				// 2^100 - 2^0	
  curve25519_sqrmodp(r4, r3);					// 2^101 - 2^1	
  curve25519_sqrmodp(r5, r4);					// 2^102 - 2^2	
  for (i = 1; i < 50; i++) {
	curve25519_sqrmodp(r4, r5);
	curve25519_sqrmodp(r5, r4);
  } 											// 2^200 - 2^100 
  curve25519_mulmodp(r4, r5, r3);				// 2^200 - 2^0	
  for (i = 0; i < 25; i++) {
	curve25519_sqrmodp(r5, r4);
	curve25519_sqrmodp(r4, r5);
  } 											// 2^250 - 2^50	
  curve25519_mulmodp(r3, r4, r2);				// 2^250 - 2^0	
  curve25519_sqrmodp(r2, r3);					// 2^251 - 2^1	
  curve25519_sqrmodp(r3, r2);					// 2^252 - 2^2	
  if (sqrtassist == 1) {
	curve25519_mulmodp(r, a, r3);				// 2^252 - 3 
  } else {
	curve25519_sqrmodp(r2, r3);					// 2^253 - 2^3	
	curve25519_sqrmodp(r3, r2);					// 2^254 - 2^4	
	curve25519_sqrmodp(r2, r3);					// 2^255 - 2^5	
	curve25519_mulmodp(r, r2, r1);				// 2^255 - 21	
  }
}

/******
 * BloodyRookie: Finding a square root mod p of x if we already know it exists and p congruent 3 mod 8.
 * Using x^((p-1)/2) congruent 1 mod p and 2^((p-1)/2) congruent -1 mod p
 * because of Eulers criterium we see that when we set v=(2x)^((p-5)/8) then
 * i:=2xv^2 is a square root of -1 and thus r=+xv(i-1) and r=-xv(i-1) are the square roots of x.
 */
curve25519_sqrtmodp = function(r, x) {
  var r1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r3 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r4 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  curve25519_addmodp(r1, x,x);								// r1 = 2x
  curve25519_invmodp(r2, r1, 1);							// r2 = (2x)^((p-5)/8) =: v
  curve25519_sqrmodp(r3, r2);								// r3 = v^2
  curve25519_mulmodp(r4, r1, r3);							// r4 = 2xv^2 =: i
  curve25519_submodp(r, r4, curve25519_one());				//  r = i-1
  curve25519_mulmodp(r1, r2, r);							// r1 = v(i-1)
  curve25519_mulmodp(r, x, r1);								//  r = xv(i-1)
}

curve25519_reduce = function (a) {
  curve25519_reduce2(a);
  
  /**
   * BloodyRookie: special case for p <= a < 2^255
   */
  if ((a[15] != 0x7FFF || a[14] != 0xFFFF || a[13] != 0xFFFF || a[12] != 0xFFFF || a[11] != 0xFFFF || a[10] != 0xFFFF || a[9] != 0xFFFF ||  a[8] != 0xFFFF || 
	   a[7] != 0xFFFF  || a[6] != 0xFFFF  || a[5] != 0xFFFF  || a[4] != 0xFFFF  || a[3] != 0xFFFF  || a[2] != 0xFFFF || a[1] != 0xFFFF || a[0] < 0xFFED))
  {
	  return;
  }
  
  var i;
  for (i=1; i<16; i++)
  {
	  a[i] = 0;
  }
  a[0] = a[0] - 0xFFED;
}
curve25519_reduce2 = function (a) {
  var v = a[15];
  if (v < 0x8000) return;
  a[15] = v % 0x8000;
  v = ~~(v / 0x8000) * 19;
  a[0] = (v += a[0]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[1] = (v += a[1]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[2] = (v += a[2]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[3] = (v += a[3]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[4] = (v += a[4]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[5] = (v += a[5]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[6] = (v += a[6]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[7] = (v += a[7]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[8] = (v += a[8]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[9] = (v += a[9]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[10] = (v += a[10]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[11] = (v += a[11]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[12] = (v += a[12]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[13] = (v += a[13]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[14] = (v += a[14]) & 0xffff;
  if ((v = ~~(v / 0x10000)) < 1) return;
  a[15] += v;
}

/**
 * Montgomery curve with A=486662 and B=1
 */
curve25519_x_to_y2 = function(r, x) {
  var r1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  curve25519_sqrmodp(r1, x);									// r1 = x^2
  curve25519_mulasmall(r2, x, 486662);							// r2 = Ax
  curve25519_addmodp(r, r1, r2);								//  r = x^2 + Ax
  curve25519_addmodp(r1, r, curve25519_one());					// r1 = x^2 + Ax + 1
  curve25519_mulmodp(r, r1, x);									//  r = x^3 + Ax^2 + x
}

curve25519_prep = function(r, s, a, b) {
  curve25519_addmodp(r, a, b);
  curve25519_submodp(s, a, b);
}

/****
 * BloodyRookie: Doubling a point on a Montgomery curve:
 * Point is given in projective coordinates p=x/z
 * 2*P = r/s, 
 * r = (x+z)^2 * (x-z)^2
 * s = ((((x+z)^2 - (x-z)^2) * 121665) + (x+z)^2) * ((x+z)^2 - (x-z)^2) 
 *   = 4*x*z * (x^2 + 486662*x*z + z^2)
 *   = 4*x*z * ((x-z)^2 + ((486662+2)/4)(4*x*z))
 */
curve25519_dbl = function(r, s, t1, t2) {
  var r1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r3 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r4 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  curve25519_sqrmodp(r1, t1);									// r1 = t1^2
  curve25519_sqrmodp(r2, t2);									// r2 = t2^2
  curve25519_submodp(r3, r1, r2);								// r3 = t1^2 - t2^2
  curve25519_mulmodp(r, r2, r1);								//  r = t1^2 * t2^2
  curve25519_mulasmall(r2, r3, 121665);							// r2 = (t1^2 - t2^2) * 121665
  curve25519_addmodp(r4, r2, r1)								// r4 = (t1^2 - t2^2) * 121665 + t1^2
  curve25519_mulmodp(s, r4, r3);								//  s = ((t1^2 - t2^2) * 121665 + t1^2) * (t1^2 - t2^2)
}

/****
 * BloodyRookie: Adding 2 points on a Montgomery curve:
 * R = Q + P = r/s when given
 * Q = x/z, P = x_p/z_p, P-Q = x_1/1
 * r = ((x-z)*(x_p+z_p) + (x+z)*(x_p-z_p))^2
 * s = x_1*((x-z)*(x_p+z_p) - (x+z)*(x_p-z_p))^2
 */
function curve25519_sum(r, s, t1, t2, t3, t4, x_1) {
  var r1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r3 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  var r4 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  curve25519_mulmodp(r1, t2, t3);								// r1 = t2 * t3
  curve25519_mulmodp(r2, t1, t4);								// r2 = t1 * t4
  curve25519_addmodp(r3, r1, r2);								// r3 = t2 * t3 + t1 * t4
  curve25519_submodp(r4, r1, r2);								// r4 = t2 * t3 - t1 * t4
  curve25519_sqrmodp(r, r3);									//  r = (t2 * t3 + t1 * t4)^2
  curve25519_sqrmodp(r1, r4);									// r1 = (t2 * t3 - t1 * t4)^2
  curve25519_mulmodp(s, r1, x_1);								//  s = (t2 * t3 - t1 * t4)^2 * x_1
}

function curve25519_(f, c, s) {
  var j, a, x_1, q, fb, counter=0; 
  var t = new Array(16), t1 = new Array(16), t2 = new Array(16), t3 = new Array(16), t4 = new Array(16);
  var sb = new Int8Array(32);
  var temp1 = new Int8Array(32);
  var temp2 = new Int8Array(64);
  var temp3 = new Int8Array(64);

  x_1 = c;
  q = [ curve25519_one(), curve25519_zero() ];
  a = [ x_1, curve25519_one() ];

  var n = 255;

  /**********************************************************************
   * BloodyRookie:                                                      *
   * Given f = f0*2^0 + f1*2^1 + ... + f255*2^255 and Basepoint a=9/1   * 
   * calculate f*a by applying the Montgomery ladder (const time algo): *
   * r0 := 0 (point at infinity)                                        *
   * r1 := a                                                            *
   * for i from 255 to 0 do                                             *
   *   if fi = 0 then                                                   *
   *      r1 := r0 + r1                                                 *          
   *      r0 := 2r0                                                     *
   *   else                                                             *
   *      r0 := r0 + r1                                                 *
   *      r1 := 2r1                                                     *
   *                                                                    *
   * Result: r0 = x-coordinate of f*a                                   *
   **********************************************************************/
  var r0 = new Array(new Array(16), new Array(16));
  var r1 = new Array(new Array(16), new Array(16));
  var t1 = new Array(16), t2 = new Array(16);
  var t3 = new Array(16), t4 = new Array(16);
  var fi;
  while (n >= 0) 
  {
    fi = curve25519_getbit(f, n);
    if (fi == 0) 
    {
       curve25519_prep(t1, t2, a[0], a[1]); 
       curve25519_prep(t3, t4, q[0], q[1]); 
       curve25519_sum(r1[0], r1[1], t1, t2, t3, t4, x_1);
       curve25519_dbl(r0[0], r0[1], t3, t4);
    } 
    else 
    {
       curve25519_prep(t1, t2, q[0], q[1]); 
       curve25519_prep(t3, t4, a[0], a[1]); 
       curve25519_sum(r0[0], r0[1], t1, t2, t3, t4, x_1);
       curve25519_dbl(r1[0], r1[1], t3, t4);
    }
    q = r0; a = r1;
    n--;
  }
  curve25519_invmodp(t, q[1], 0);
  curve25519_mulmodp(t1, q[0], t);
  q[0] = curve25519_cpy16(t1);
  
  // q[0]=x-coordinate of k*G=:Px
  // q[1]=z-coordinate of k*G=:Pz
  // a = q + G = P + G
  if (s != null)
  {
	  /*************************************************************************
	   * BloodyRookie: Recovery of the y-coordinate of point P:                *
	   *                                                                       *
	   * If P=(x,y), P1=(x1, y1), P2=(x2,y2) and P2 = P1 + P then              *
	   *                                                                       *
	   * y1 = ((x1 * x + 1)(x1 + x + 2A) - 2A - (x1 - x)^2 * x2)/2y            *
	   *                                                                       *
	   * Setting P2=Q, P1=P and P=G in the above formula we get                *
	   *                                                                       *
	   * Py =  ((Px * Gx + 1) * (Px + Gx + 2A) - 2A - (Px - Gx)^2 * Qx)/(2*Gy) *
	   *    = -((Qx + Px + Gx + A) * (Px - Gx)^2 - Py^2 - Gy^2)/(2*Gy)         *
	   *************************************************************************/
	  t = curve25519_cpy16(q[0]);
	  curve25519_x_to_y2(t1, t);								// t1 = Py^2
	  curve25519_invmodp(t3, a[1], 0);
	  curve25519_mulmodp(t2, a[0], t3);							// t2 = (P+G)x = Qx
	  curve25519_addmodp(t4, t2, t);							// t4 =  Qx + Px
	  curve25519_addmodp(t2, t4, curve25519_486671());			// t2 = Qx + Px + Gx + A  
	  curve25519_submodp(t4, t, curve25519_nine());				// t4 = Px - Gx
	  curve25519_sqrmodp(t3, t4);								// t3 = (Px - Gx)^2
	  curve25519_mulmodp(t4, t2, t3);							// t4 = (Qx + Px + Gx + A) * (Px - Gx)^2
	  curve25519_submodp(t, t4, t1);							//  t = (Qx + Px + Gx + A) * (Px - Gx)^2 - Py^2
	  curve25519_submodp(t4, t, curve25519_39420360());			// t4 = (Qx + Px + Gx + A) * (Px - Gx)^2 - Py^2 - Gy^2
	  curve25519_mulmodp(t1, t4, curve25519_r2y())				// t1 = ((Qx + Px + Gx + A) * (Px - Gx)^2 - Py^2 - Gy^2)/(2Gy) = -Py
	  fb = curve25519_convertToByteArray(f);
	  j = curve25519_isNegative(t1);
	  if (j != 0)
	  {
		  /***
		   * Py is positiv, so just copy
		   */
		  sb = curve25519_cpy32(fb);
	  }
	  else
	  {
		  /***
		   * Py is negative:
		   * We will take s = -f^-1 mod q instead of s=f^-1 mod q
		   */
		  curve25519_mula_small(sb, curve25519_order_times_8, 0, fb, 32, -1);
	  }
	  
	  temp1 = curve25519_cpy32(curve25519_order);
	  temp1 = curve25519_egcd32(temp2, temp3, sb, temp1);
	  sb = curve25519_cpy32(temp1);
	  if ((sb[31] & 0x80)!=0)
	  {
		  curve25519_mula_small(sb, sb, 0, curve25519_order, 32, 1);
	  }
	  var stmp = curve25519_convertToShortArray(sb);
	  curve25519_fillShortArray(stmp, s);
  }

  return q[0];
}

curve25519_keygen = function(s, curve) {
	curve25519_clamp(curve);
	return curve25519_(curve, curve25519_nine(), s);
}

/* Signature generation primitive, calculates (x-h)s mod q
 *   v  [out] signature value
 *   h  [in]  signature hash (of message, signature pub key, and context data)
 *   x  [in]  signature private key
 *   s  [in]  private key for signing
 * returns true on success, false on failure (use different x or h)
 */
curve25519_sign = function(v, h, x, s) {
	tmp1=new Int8Array(65);
	tmp2=new Int8Array(33);
	for (i = 0; i < 32; i++)
	{
		v[i] = 0;
	}
	curve25519_mula_small(v, x, 0, h, 32, -1);
	curve25519_mula_small(v, v, 0, curve25519_order, 32, parseInt((15-v[31])/16));
	curve25519_mula32(tmp1, v, s, 32, 1);
	curve25519_divmod(tmp2, tmp1, 64, curve25519_order, 32);
	w=0;
	for (k = 0; k < 32; k++)
	{
		v[k] = tmp1[k];
		w |= v[k];
	}
	return w != 0;
}

curve25519_verify = function(Y, v, h, P) {
	d=new Int8Array(32);
	yx=new Array(new Array(16), new Array(16), new Array(16));
	yz=new Array(new Array(16), new Array(16), new Array(16));
	var s=new Array(new Array(16), new Array(16));
	var q=new Array(new Array(16), new Array(16));
	var t1=new Array(new Array(16), new Array(16), new Array(16));
	var t2=new Array(new Array(16), new Array(16), new Array(16));
	var vi = 0, hi = 0, di = 0, nvh=0, i=0, j=0, k=0, counter=1;

	/******************************************************************
     * Set s[0] to P+G and s[1] to P-G.                               *
     * If sqrt(Py^2) is negativ we switch s[0] and s[1]               *
	 *                                                                *
     * s[0] = (Py^2 + Gy^2 - 2 Py Gy)/(Px - Gx)^2 - Px - Gx - 486662  *
     * s[1] = (Py^2 + Gy^2 + 2 Py Gy)/(Px - Gx)^2 - Px - Gx - 486662  *
     ******************************************************************/

	var p = [ curve25519_nine(), curve25519_convertToShortArray(P) ];
	curve25519_x_to_y2(q[0], p[1]);								// q[0] = Py^2
	curve25519_sqrtmodp(t1[0], q[0]);							// t1[0] = +-Py
	j = curve25519_isNegative(t1[0]);
	curve25519_addmodp(t2[0], q[0], curve25519_39420360());		// t2[0] = Py^2 + Gy^2
	curve25519_mulmodp(t2[1], curve25519_2y(), t1[0]);			// t2[1] = +-Py * 2Gy
	curve25519_submodp(t1[j], t2[0], t2[1]);					// t1[j] = Py^2 + Gy^2 - +-Py * 2Gy
	curve25519_addmodp(t1[1-j], t2[0], t2[1]);					// t1[1-j] = Py^2 + Gy^2 + +-Py * 2Gy
	q[0] = curve25519_cpy16(p[1]);								// q[0] = Px
	curve25519_submodp(t2[0], q[0], curve25519_nine());			// t2[0] = Px-Gx
	curve25519_sqrmodp(t2[1], t2[0]);							// t2[1] = (Px-Gx)^2
	curve25519_invmodp(t2[0], t2[1], 0);						// t2[0] = 1/(Px-Gx)^2
	curve25519_mulmodp(q[0], t1[0], t2[0]);						// q[0] = (Py^2 + Gy^2 - Py * 2Gy)/(Px-Gx)^2
	curve25519_submodp(q[1], q[0], p[1]);						// q[1] = (Py^2 + Gy^2 - Py * 2Gy)/(Px-Gx)^2 - Px
	curve25519_submodp(s[0], q[1], curve25519_486671());		// s[0] = (Py^2 + Gy^2 - Py * 2Gy)/(Px-Gx)^2 - Px - Gx - A = P+Q
	curve25519_mulmodp(q[0], t1[1], t2[0]);						// q[0] = (Py^2 + Gy^2 + Py * 2Gy)/(Px-Gx)^2
	curve25519_submodp(q[1], q[0], p[1]);						// q[1] = (Py^2 + Gy^2 + Py * 2Gy)/(Px-Gx)^2 - Px
	curve25519_submodp(s[1], q[1], curve25519_486671());		// s[1] = (Py^2 + Gy^2 + Py * 2Gy)/(Px-Gx)^2 - Px - Gx - A = P-Q
	
	/**
	 * Fast algorithm for computing vP+hG
	 */
	for (i = 0; i < 32; i++) 
	{
		vi = (vi >> 8) ^ (v[i] & 0xFF) ^ ((v[i] & 0xFF) << 1);
		hi = (hi >> 8) ^ (h[i] & 0xFF) ^ ((h[i] & 0xFF) << 1);
		nvh = ~(vi ^ hi);
		di = (nvh & (di & 0x80) >> 7) ^ vi;
		di ^= nvh & (di & 0x01) << 1;
		di ^= nvh & (di & 0x02) << 1;
		di ^= nvh & (di & 0x04) << 1;
		di ^= nvh & (di & 0x08) << 1;
		di ^= nvh & (di & 0x10) << 1;
		di ^= nvh & (di & 0x20) << 1;
		di ^= nvh & (di & 0x40) << 1;
		d[i] = (di & 0xFF);
	}

	di = ((nvh & (di & 0x80) << 1) ^ vi) >> 8;

	/**
	 * yx[0]/yz[0] = point at infinity
	 */
	yx[0] = curve25519_cpy16(curve25519_one());
	yx[1] = curve25519_cpy16(p[di]);
	yx[2] = curve25519_cpy16(s[0]);
	yz[0] = curve25519_cpy16(curve25519_zero());
	yz[1] = curve25519_cpy16(curve25519_one());
	yz[2] = curve25519_cpy16(curve25519_one());
	
	vi = 0;
	hi = 0;

	for (i = 32; i-- != 0; i=i) 
	{
		vi = (vi << 8) | (v[i] & 0xFF);
		hi = (hi << 8) | (h[i] & 0xFF);
		di = (di << 8) | (d[i] & 0xFF);

		for (j = 8; j-- !=0 ; j=j) 
		{
			k = ((vi ^ vi >> 1) >> j & 1) + ((hi ^ hi >> 1) >> j & 1);
			curve25519_prep(t1[0], t2[0], yx[0], yz[0]);
			curve25519_prep(t1[1], t2[1], yx[1], yz[1]);
			curve25519_prep(t1[2], t2[2], yx[2], yz[2]);
			
			curve25519_dbl(yx[0], yz[0], t1[k], t2[k]);
			k = (di >> j & 2) ^ ((di >> j & 1) << 1);
			curve25519_sum(yx[1], yz[1], t1[1], t2[1], t1[k], t2[k], p[di >> j & 1]);
			curve25519_sum(yx[2], yz[2], t1[2], t2[2], t1[0], t2[0], s[((vi ^ hi) >> j & 2) >> 1]);
		}
	}

	k = (vi & 1) + (hi & 1);
	curve25519_invmodp(t1[0], yz[k], 0);
	curve25519_mulmodp(t1[1], yx[k], t1[0]);
	var YY = curve25519_convertToByteArray(t1[1]);
	curve25519_fillByteArray(YY, Y);
}

// Copyright (c) 2005  Tom Wu
// All Rights Reserved.
// See "LICENSE" for details.

// Basic JavaScript BN library - subset useful for RSA encryption.

// Bits per digit
var dbits;

// JavaScript engine analysis
var canary = 0xdeadbeefcafe;
var j_lm = ((canary&0xffffff)==0xefcafe);

// (public) Constructor
function BigInteger(a,b,c) {
  if(a != null)
    if("number" == typeof a) this.fromNumber(a,b,c);
    else if(b == null && "string" != typeof a) this.fromString(a,256);
    else this.fromString(a,b);
}

// return new, unset BigInteger
function nbi() { return new BigInteger(null); }

// am: Compute w_j += (x*this_i), propagate carries,
// c is initial carry, returns final carry.
// c < 3*dvalue, x < 2*dvalue, this_i < dvalue
// We need to select the fastest one that works in this environment.

// am1: use a single mult and divide to get the high bits,
// max digit bits should be 26 because
// max internal value = 2*dvalue^2-2*dvalue (< 2^53)
function am1(i,x,w,j,c,n) {
  while(--n >= 0) {
    var v = x*this[i++]+w[j]+c;
    c = Math.floor(v/0x4000000);
    w[j++] = v&0x3ffffff;
  }
  return c;
}
// am2 avoids a big mult-and-extract completely.
// Max digit bits should be <= 30 because we do bitwise ops
// on values up to 2*hdvalue^2-hdvalue-1 (< 2^31)
function am2(i,x,w,j,c,n) {
  var xl = x&0x7fff, xh = x>>15;
  while(--n >= 0) {
    var l = this[i]&0x7fff;
    var h = this[i++]>>15;
    var m = xh*l+h*xl;
    l = xl*l+((m&0x7fff)<<15)+w[j]+(c&0x3fffffff);
    c = (l>>>30)+(m>>>15)+xh*h+(c>>>30);
    w[j++] = l&0x3fffffff;
  }
  return c;
}
// Alternately, set max digit bits to 28 since some
// browsers slow down when dealing with 32-bit numbers.
function am3(i,x,w,j,c,n) {
  var xl = x&0x3fff, xh = x>>14;
  while(--n >= 0) {
    var l = this[i]&0x3fff;
    var h = this[i++]>>14;
    var m = xh*l+h*xl;
    l = xl*l+((m&0x3fff)<<14)+w[j]+c;
    c = (l>>28)+(m>>14)+xh*h;
    w[j++] = l&0xfffffff;
  }
  return c;
}

q = [ appName = "Microsoft Internet Explorer" ]

if(j_lm && (q.appName == "Microsoft Internet Explorer")) {
  BigInteger.prototype.am = am2;
  dbits = 30;
}
else if(j_lm && (q.appName != "Netscape")) {
  BigInteger.prototype.am = am1;
  dbits = 26;
}
else { // Mozilla/Netscape seems to prefer am3
  BigInteger.prototype.am = am3;
  dbits = 28;
}

BigInteger.prototype.DB = dbits;
BigInteger.prototype.DM = ((1<<dbits)-1);
BigInteger.prototype.DV = (1<<dbits);

var BI_FP = 52;
BigInteger.prototype.FV = Math.pow(2,BI_FP);
BigInteger.prototype.F1 = BI_FP-dbits;
BigInteger.prototype.F2 = 2*dbits-BI_FP;

// Digit conversions
var BI_RM = "0123456789abcdefghijklmnopqrstuvwxyz";
var BI_RC = new Array();
var rr,vv;
rr = "0".charCodeAt(0);
for(vv = 0; vv <= 9; ++vv) BI_RC[rr++] = vv;
rr = "a".charCodeAt(0);
for(vv = 10; vv < 36; ++vv) BI_RC[rr++] = vv;
rr = "A".charCodeAt(0);
for(vv = 10; vv < 36; ++vv) BI_RC[rr++] = vv;

function int2char(n) { return BI_RM.charAt(n); }
function intAt(s,i) {
  var c = BI_RC[s.charCodeAt(i)];
  return (c==null)?-1:c;
}

// (protected) copy this to r
function bnpCopyTo(r) {
  for(var i = this.t-1; i >= 0; --i) r[i] = this[i];
  r.t = this.t;
  r.s = this.s;
}

// (protected) set from integer value x, -DV <= x < DV
function bnpFromInt(x) {
  this.t = 1;
  this.s = (x<0)?-1:0;
  if(x > 0) this[0] = x;
  else if(x < -1) this[0] = x+this.DV;
  else this.t = 0;
}

// return bigint initialized to value
function nbv(i) { var r = nbi(); r.fromInt(i); return r; }

// (protected) set from string and radix
function bnpFromString(s,b) {
  var k;
  if(b == 16) k = 4;
  else if(b == 8) k = 3;
  else if(b == 256) k = 8; // byte array
  else if(b == 2) k = 1;
  else if(b == 32) k = 5;
  else if(b == 4) k = 2;
  else { this.fromRadix(s,b); return; }
  this.t = 0;
  this.s = 0;
  var i = s.length, mi = false, sh = 0;
  while(--i >= 0) {
    var x = (k==8)?s[i]&0xff:intAt(s,i);
    if(x < 0) {
      if(s.charAt(i) == "-") mi = true;
      continue;
    }
    mi = false;
    if(sh == 0)
      this[this.t++] = x;
    else if(sh+k > this.DB) {
      this[this.t-1] |= (x&((1<<(this.DB-sh))-1))<<sh;
      this[this.t++] = (x>>(this.DB-sh));
    }
    else
      this[this.t-1] |= x<<sh;
    sh += k;
    if(sh >= this.DB) sh -= this.DB;
  }
  if(k == 8 && (s[0]&0x80) != 0) {
    this.s = -1;
    if(sh > 0) this[this.t-1] |= ((1<<(this.DB-sh))-1)<<sh;
  }
  this.clamp();
  if(mi) BigInteger.ZERO.subTo(this,this);
}

// (protected) clamp off excess high words
function bnpClamp() {
  var c = this.s&this.DM;
  while(this.t > 0 && this[this.t-1] == c) --this.t;
}

// (public) return string representation in given radix
function bnToString(b) {
  if(this.s < 0) return "-"+this.negate().toString(b);
  var k;
  if(b == 16) k = 4;
  else if(b == 8) k = 3;
  else if(b == 2) k = 1;
  else if(b == 32) k = 5;
  else if(b == 4) k = 2;
  else return this.toRadix(b);
  var km = (1<<k)-1, d, m = false, r = "", i = this.t;
  var p = this.DB-(i*this.DB)%k;
  if(i-- > 0) {
    if(p < this.DB && (d = this[i]>>p) > 0) { m = true; r = int2char(d); }
    while(i >= 0) {
      if(p < k) {
        d = (this[i]&((1<<p)-1))<<(k-p);
        d |= this[--i]>>(p+=this.DB-k);
      }
      else {
        d = (this[i]>>(p-=k))&km;
        if(p <= 0) { p += this.DB; --i; }
      }
      if(d > 0) m = true;
      if(m) r += int2char(d);
    }
  }
  return m?r:"0";
}

// (public) -this
function bnNegate() { var r = nbi(); BigInteger.ZERO.subTo(this,r); return r; }

// (public) |this|
function bnAbs() { return (this.s<0)?this.negate():this; }

// (public) return + if this > a, - if this < a, 0 if equal
function bnCompareTo(a) {
  var r = this.s-a.s;
  if(r != 0) return r;
  var i = this.t;
  r = i-a.t;
  if(r != 0) return (this.s<0)?-r:r;
  while(--i >= 0) if((r=this[i]-a[i]) != 0) return r;
  return 0;
}

// returns bit length of the integer x
function nbits(x) {
  var r = 1, t;
  if((t=x>>>16) != 0) { x = t; r += 16; }
  if((t=x>>8) != 0) { x = t; r += 8; }
  if((t=x>>4) != 0) { x = t; r += 4; }
  if((t=x>>2) != 0) { x = t; r += 2; }
  if((t=x>>1) != 0) { x = t; r += 1; }
  return r;
}

// (public) return the number of bits in "this"
function bnBitLength() {
  if(this.t <= 0) return 0;
  return this.DB*(this.t-1)+nbits(this[this.t-1]^(this.s&this.DM));
}

// (protected) r = this << n*DB
function bnpDLShiftTo(n,r) {
  var i;
  for(i = this.t-1; i >= 0; --i) r[i+n] = this[i];
  for(i = n-1; i >= 0; --i) r[i] = 0;
  r.t = this.t+n;
  r.s = this.s;
}

// (protected) r = this >> n*DB
function bnpDRShiftTo(n,r) {
  for(var i = n; i < this.t; ++i) r[i-n] = this[i];
  r.t = Math.max(this.t-n,0);
  r.s = this.s;
}

// (protected) r = this << n
function bnpLShiftTo(n,r) {
  var bs = n%this.DB;
  var cbs = this.DB-bs;
  var bm = (1<<cbs)-1;
  var ds = Math.floor(n/this.DB), c = (this.s<<bs)&this.DM, i;
  for(i = this.t-1; i >= 0; --i) {
    r[i+ds+1] = (this[i]>>cbs)|c;
    c = (this[i]&bm)<<bs;
  }
  for(i = ds-1; i >= 0; --i) r[i] = 0;
  r[ds] = c;
  r.t = this.t+ds+1;
  r.s = this.s;
  r.clamp();
}

// (protected) r = this >> n
function bnpRShiftTo(n,r) {
  r.s = this.s;
  var ds = Math.floor(n/this.DB);
  if(ds >= this.t) { r.t = 0; return; }
  var bs = n%this.DB;
  var cbs = this.DB-bs;
  var bm = (1<<bs)-1;
  r[0] = this[ds]>>bs;
  for(var i = ds+1; i < this.t; ++i) {
    r[i-ds-1] |= (this[i]&bm)<<cbs;
    r[i-ds] = this[i]>>bs;
  }
  if(bs > 0) r[this.t-ds-1] |= (this.s&bm)<<cbs;
  r.t = this.t-ds;
  r.clamp();
}

// (protected) r = this - a
function bnpSubTo(a,r) {
  var i = 0, c = 0, m = Math.min(a.t,this.t);
  while(i < m) {
    c += this[i]-a[i];
    r[i++] = c&this.DM;
    c >>= this.DB;
  }
  if(a.t < this.t) {
    c -= a.s;
    while(i < this.t) {
      c += this[i];
      r[i++] = c&this.DM;
      c >>= this.DB;
    }
    c += this.s;
  }
  else {
    c += this.s;
    while(i < a.t) {
      c -= a[i];
      r[i++] = c&this.DM;
      c >>= this.DB;
    }
    c -= a.s;
  }
  r.s = (c<0)?-1:0;
  if(c < -1) r[i++] = this.DV+c;
  else if(c > 0) r[i++] = c;
  r.t = i;
  r.clamp();
}

// (protected) r = this * a, r != this,a (HAC 14.12)
// "this" should be the larger one if appropriate.
function bnpMultiplyTo(a,r) {
  var x = this.abs(), y = a.abs();
  var i = x.t;
  r.t = i+y.t;
  while(--i >= 0) r[i] = 0;
  for(i = 0; i < y.t; ++i) r[i+x.t] = x.am(0,y[i],r,i,0,x.t);
  r.s = 0;
  r.clamp();
  if(this.s != a.s) BigInteger.ZERO.subTo(r,r);
}

// (protected) r = this^2, r != this (HAC 14.16)
function bnpSquareTo(r) {
  var x = this.abs();
  var i = r.t = 2*x.t;
  while(--i >= 0) r[i] = 0;
  for(i = 0; i < x.t-1; ++i) {
    var c = x.am(i,x[i],r,2*i,0,1);
    if((r[i+x.t]+=x.am(i+1,2*x[i],r,2*i+1,c,x.t-i-1)) >= x.DV) {
      r[i+x.t] -= x.DV;
      r[i+x.t+1] = 1;
    }
  }
  if(r.t > 0) r[r.t-1] += x.am(i,x[i],r,2*i,0,1);
  r.s = 0;
  r.clamp();
}

// (protected) divide this by m, quotient and remainder to q, r (HAC 14.20)
// r != q, this != m.  q or r may be null.
function bnpDivRemTo(m,q,r) {
  var pm = m.abs();
  if(pm.t <= 0) return;
  var pt = this.abs();
  if(pt.t < pm.t) {
    if(q != null) q.fromInt(0);
    if(r != null) this.copyTo(r);
    return;
  }
  if(r == null) r = nbi();
  var y = nbi(), ts = this.s, ms = m.s;
  var nsh = this.DB-nbits(pm[pm.t-1]);	// normalize modulus
  if(nsh > 0) { pm.lShiftTo(nsh,y); pt.lShiftTo(nsh,r); }
  else { pm.copyTo(y); pt.copyTo(r); }
  var ys = y.t;
  var y0 = y[ys-1];
  if(y0 == 0) return;
  var yt = y0*(1<<this.F1)+((ys>1)?y[ys-2]>>this.F2:0);
  var d1 = this.FV/yt, d2 = (1<<this.F1)/yt, e = 1<<this.F2;
  var i = r.t, j = i-ys, t = (q==null)?nbi():q;
  y.dlShiftTo(j,t);
  if(r.compareTo(t) >= 0) {
    r[r.t++] = 1;
    r.subTo(t,r);
  }
  BigInteger.ONE.dlShiftTo(ys,t);
  t.subTo(y,y);	// "negative" y so we can replace sub with am later
  while(y.t < ys) y[y.t++] = 0;
  while(--j >= 0) {
    // Estimate quotient digit
    var qd = (r[--i]==y0)?this.DM:Math.floor(r[i]*d1+(r[i-1]+e)*d2);
    if((r[i]+=y.am(0,qd,r,j,0,ys)) < qd) {	// Try it out
      y.dlShiftTo(j,t);
      r.subTo(t,r);
      while(r[i] < --qd) r.subTo(t,r);
    }
  }
  if(q != null) {
    r.drShiftTo(ys,q);
    if(ts != ms) BigInteger.ZERO.subTo(q,q);
  }
  r.t = ys;
  r.clamp();
  if(nsh > 0) r.rShiftTo(nsh,r);	// Denormalize remainder
  if(ts < 0) BigInteger.ZERO.subTo(r,r);
}

// (public) this mod a
function bnMod(a) {
  var r = nbi();
  this.abs().divRemTo(a,null,r);
  if(this.s < 0 && r.compareTo(BigInteger.ZERO) > 0) a.subTo(r,r);
  return r;
}

// Modular reduction using "classic" algorithm
function Classic(m) { this.m = m; }
function cConvert(x) {
  if(x.s < 0 || x.compareTo(this.m) >= 0) return x.mod(this.m);
  else return x;
}
function cRevert(x) { return x; }
function cReduce(x) { x.divRemTo(this.m,null,x); }
function cMulTo(x,y,r) { x.multiplyTo(y,r); this.reduce(r); }
function cSqrTo(x,r) { x.squareTo(r); this.reduce(r); }

Classic.prototype.convert = cConvert;
Classic.prototype.revert = cRevert;
Classic.prototype.reduce = cReduce;
Classic.prototype.mulTo = cMulTo;
Classic.prototype.sqrTo = cSqrTo;

// (protected) return "-1/this % 2^DB"; useful for Mont. reduction
// justification:
//         xy == 1 (mod m)
//         xy =  1+km
//   xy(2-xy) = (1+km)(1-km)
// x[y(2-xy)] = 1-k^2m^2
// x[y(2-xy)] == 1 (mod m^2)
// if y is 1/x mod m, then y(2-xy) is 1/x mod m^2
// should reduce x and y(2-xy) by m^2 at each step to keep size bounded.
// JS multiply "overflows" differently from C/C++, so care is needed here.
function bnpInvDigit() {
  if(this.t < 1) return 0;
  var x = this[0];
  if((x&1) == 0) return 0;
  var y = x&3;		// y == 1/x mod 2^2
  y = (y*(2-(x&0xf)*y))&0xf;	// y == 1/x mod 2^4
  y = (y*(2-(x&0xff)*y))&0xff;	// y == 1/x mod 2^8
  y = (y*(2-(((x&0xffff)*y)&0xffff)))&0xffff;	// y == 1/x mod 2^16
  // last step - calculate inverse mod DV directly;
  // assumes 16 < DB <= 32 and assumes ability to handle 48-bit ints
  y = (y*(2-x*y%this.DV))%this.DV;		// y == 1/x mod 2^dbits
  // we really want the negative inverse, and -DV < y < DV
  return (y>0)?this.DV-y:-y;
}

// Montgomery reduction
function Montgomery(m) {
  this.m = m;
  this.mp = m.invDigit();
  this.mpl = this.mp&0x7fff;
  this.mph = this.mp>>15;
  this.um = (1<<(m.DB-15))-1;
  this.mt2 = 2*m.t;
}

// xR mod m
function montConvert(x) {
  var r = nbi();
  x.abs().dlShiftTo(this.m.t,r);
  r.divRemTo(this.m,null,r);
  if(x.s < 0 && r.compareTo(BigInteger.ZERO) > 0) this.m.subTo(r,r);
  return r;
}

// x/R mod m
function montRevert(x) {
  var r = nbi();
  x.copyTo(r);
  this.reduce(r);
  return r;
}

// x = x/R mod m (HAC 14.32)
function montReduce(x) {
  while(x.t <= this.mt2)	// pad x so am has enough room later
    x[x.t++] = 0;
  for(var i = 0; i < this.m.t; ++i) {
    // faster way of calculating u0 = x[i]*mp mod DV
    var j = x[i]&0x7fff;
    var u0 = (j*this.mpl+(((j*this.mph+(x[i]>>15)*this.mpl)&this.um)<<15))&x.DM;
    // use am to combine the multiply-shift-add into one call
    j = i+this.m.t;
    x[j] += this.m.am(0,u0,x,i,0,this.m.t);
    // propagate carry
    while(x[j] >= x.DV) { x[j] -= x.DV; x[++j]++; }
  }
  x.clamp();
  x.drShiftTo(this.m.t,x);
  if(x.compareTo(this.m) >= 0) x.subTo(this.m,x);
}

// r = "x^2/R mod m"; x != r
function montSqrTo(x,r) { x.squareTo(r); this.reduce(r); }

// r = "xy/R mod m"; x,y != r
function montMulTo(x,y,r) { x.multiplyTo(y,r); this.reduce(r); }

Montgomery.prototype.convert = montConvert;
Montgomery.prototype.revert = montRevert;
Montgomery.prototype.reduce = montReduce;
Montgomery.prototype.mulTo = montMulTo;
Montgomery.prototype.sqrTo = montSqrTo;

// (protected) true iff this is even
function bnpIsEven() { return ((this.t>0)?(this[0]&1):this.s) == 0; }

// (protected) this^e, e < 2^32, doing sqr and mul with "r" (HAC 14.79)
function bnpExp(e,z) {
  if(e > 0xffffffff || e < 1) return BigInteger.ONE;
  var r = nbi(), r2 = nbi(), g = z.convert(this), i = nbits(e)-1;
  g.copyTo(r);
  while(--i >= 0) {
    z.sqrTo(r,r2);
    if((e&(1<<i)) > 0) z.mulTo(r2,g,r);
    else { var t = r; r = r2; r2 = t; }
  }
  return z.revert(r);
}

// (public) this^e % m, 0 <= e < 2^32
function bnModPowInt(e,m) {
  var z;
  if(e < 256 || m.isEven()) z = new Classic(m); else z = new Montgomery(m);
  return this.exp(e,z);
}

// protected
BigInteger.prototype.copyTo = bnpCopyTo;
BigInteger.prototype.fromInt = bnpFromInt;
BigInteger.prototype.fromString = bnpFromString;
BigInteger.prototype.clamp = bnpClamp;
BigInteger.prototype.dlShiftTo = bnpDLShiftTo;
BigInteger.prototype.drShiftTo = bnpDRShiftTo;
BigInteger.prototype.lShiftTo = bnpLShiftTo;
BigInteger.prototype.rShiftTo = bnpRShiftTo;
BigInteger.prototype.subTo = bnpSubTo;
BigInteger.prototype.multiplyTo = bnpMultiplyTo;
BigInteger.prototype.squareTo = bnpSquareTo;
BigInteger.prototype.divRemTo = bnpDivRemTo;
BigInteger.prototype.invDigit = bnpInvDigit;
BigInteger.prototype.isEven = bnpIsEven;
BigInteger.prototype.exp = bnpExp;

// public
BigInteger.prototype.toString = bnToString;
BigInteger.prototype.negate = bnNegate;
BigInteger.prototype.abs = bnAbs;
BigInteger.prototype.compareTo = bnCompareTo;
BigInteger.prototype.bitLength = bnBitLength;
BigInteger.prototype.mod = bnMod;
BigInteger.prototype.modPowInt = bnModPowInt;

// "constants"
BigInteger.ZERO = nbv(0);
BigInteger.ONE = nbv(1);

// Copyright (c) 2005-2009  Tom Wu
// All Rights Reserved.
// See "LICENSE" for details.

// Extended JavaScript BN functions, required for RSA private ops.

// Version 1.1: new BigInteger("0", 10) returns "proper" zero
// Version 1.2: square() API, isProbablePrime fix

// (public)
function bnClone() { var r = nbi(); this.copyTo(r); return r; }

// (public) return value as integer
function bnIntValue() {
  if(this.s < 0) {
    if(this.t == 1) return this[0]-this.DV;
    else if(this.t == 0) return -1;
  }
  else if(this.t == 1) return this[0];
  else if(this.t == 0) return 0;
  // assumes 16 < DB < 32
  return ((this[1]&((1<<(32-this.DB))-1))<<this.DB)|this[0];
}

// (public) return value as byte
function bnByteValue() { return (this.t==0)?this.s:(this[0]<<24)>>24; }

// (public) return value as short (assumes DB>=16)
function bnShortValue() { return (this.t==0)?this.s:(this[0]<<16)>>16; }

// (protected) return x s.t. r^x < DV
function bnpChunkSize(r) { return Math.floor(Math.LN2*this.DB/Math.log(r)); }

// (public) 0 if this == 0, 1 if this > 0
function bnSigNum() {
  if(this.s < 0) return -1;
  else if(this.t <= 0 || (this.t == 1 && this[0] <= 0)) return 0;
  else return 1;
}

// (protected) convert to radix string
function bnpToRadix(b) {
  if(b == null) b = 10;
  if(this.signum() == 0 || b < 2 || b > 36) return "0";
  var cs = this.chunkSize(b);
  var a = Math.pow(b,cs);
  var d = nbv(a), y = nbi(), z = nbi(), r = "";
  this.divRemTo(d,y,z);
  while(y.signum() > 0) {
    r = (a+z.intValue()).toString(b).substr(1) + r;
    y.divRemTo(d,y,z);
  }
  return z.intValue().toString(b) + r;
}

// (protected) convert from radix string
function bnpFromRadix(s,b) {
  this.fromInt(0);
  if(b == null) b = 10;
  var cs = this.chunkSize(b);
  var d = Math.pow(b,cs), mi = false, j = 0, w = 0;
  for(var i = 0; i < s.length; ++i) {
    var x = intAt(s,i);
    if(x < 0) {
      if(s.charAt(i) == "-" && this.signum() == 0) mi = true;
      continue;
    }
    w = b*w+x;
    if(++j >= cs) {
      this.dMultiply(d);
      this.dAddOffset(w,0);
      j = 0;
      w = 0;
    }
  }
  if(j > 0) {
    this.dMultiply(Math.pow(b,j));
    this.dAddOffset(w,0);
  }
  if(mi) BigInteger.ZERO.subTo(this,this);
}

// (protected) alternate constructor
function bnpFromNumber(a,b,c) {
  if("number" == typeof b) {
    // new BigInteger(int,int,RNG)
    if(a < 2) this.fromInt(1);
    else {
      this.fromNumber(a,c);
      if(!this.testBit(a-1))	// force MSB set
        this.bitwiseTo(BigInteger.ONE.shiftLeft(a-1),op_or,this);
      if(this.isEven()) this.dAddOffset(1,0); // force odd
      while(!this.isProbablePrime(b)) {
        this.dAddOffset(2,0);
        if(this.bitLength() > a) this.subTo(BigInteger.ONE.shiftLeft(a-1),this);
      }
    }
  }
  else {
    // new BigInteger(int,RNG)
    var x = new Array(), t = a&7;
    x.length = (a>>3)+1;
    b.nextBytes(x);
    if(t > 0) x[0] &= ((1<<t)-1); else x[0] = 0;
    this.fromString(x,256);
  }
}

// (public) convert to bigendian byte array
function bnToByteArray() {
  var i = this.t, r = new Array();
  r[0] = this.s;
  var p = this.DB-(i*this.DB)%8, d, k = 0;
  if(i-- > 0) {
    if(p < this.DB && (d = this[i]>>p) != (this.s&this.DM)>>p)
      r[k++] = d|(this.s<<(this.DB-p));
    while(i >= 0) {
      if(p < 8) {
        d = (this[i]&((1<<p)-1))<<(8-p);
        d |= this[--i]>>(p+=this.DB-8);
      }
      else {
        d = (this[i]>>(p-=8))&0xff;
        if(p <= 0) { p += this.DB; --i; }
      }
      if((d&0x80) != 0) d |= -256;
      if(k == 0 && (this.s&0x80) != (d&0x80)) ++k;
      if(k > 0 || d != this.s) r[k++] = d;
    }
  }
  return r;
}

function bnEquals(a) { return(this.compareTo(a)==0); }
function bnMin(a) { return(this.compareTo(a)<0)?this:a; }
function bnMax(a) { return(this.compareTo(a)>0)?this:a; }

// (protected) r = this op a (bitwise)
function bnpBitwiseTo(a,op,r) {
  var i, f, m = Math.min(a.t,this.t);
  for(i = 0; i < m; ++i) r[i] = op(this[i],a[i]);
  if(a.t < this.t) {
    f = a.s&this.DM;
    for(i = m; i < this.t; ++i) r[i] = op(this[i],f);
    r.t = this.t;
  }
  else {
    f = this.s&this.DM;
    for(i = m; i < a.t; ++i) r[i] = op(f,a[i]);
    r.t = a.t;
  }
  r.s = op(this.s,a.s);
  r.clamp();
}

// (public) this & a
function op_and(x,y) { return x&y; }
function bnAnd(a) { var r = nbi(); this.bitwiseTo(a,op_and,r); return r; }

// (public) this | a
function op_or(x,y) { return x|y; }
function bnOr(a) { var r = nbi(); this.bitwiseTo(a,op_or,r); return r; }

// (public) this ^ a
function op_xor(x,y) { return x^y; }
function bnXor(a) { var r = nbi(); this.bitwiseTo(a,op_xor,r); return r; }

// (public) this & ~a
function op_andnot(x,y) { return x&~y; }
function bnAndNot(a) { var r = nbi(); this.bitwiseTo(a,op_andnot,r); return r; }

// (public) ~this
function bnNot() {
  var r = nbi();
  for(var i = 0; i < this.t; ++i) r[i] = this.DM&~this[i];
  r.t = this.t;
  r.s = ~this.s;
  return r;
}

// (public) this << n
function bnShiftLeft(n) {
  var r = nbi();
  if(n < 0) this.rShiftTo(-n,r); else this.lShiftTo(n,r);
  return r;
}

// (public) this >> n
function bnShiftRight(n) {
  var r = nbi();
  if(n < 0) this.lShiftTo(-n,r); else this.rShiftTo(n,r);
  return r;
}

// return index of lowest 1-bit in x, x < 2^31
function lbit(x) {
  if(x == 0) return -1;
  var r = 0;
  if((x&0xffff) == 0) { x >>= 16; r += 16; }
  if((x&0xff) == 0) { x >>= 8; r += 8; }
  if((x&0xf) == 0) { x >>= 4; r += 4; }
  if((x&3) == 0) { x >>= 2; r += 2; }
  if((x&1) == 0) ++r;
  return r;
}

// (public) returns index of lowest 1-bit (or -1 if none)
function bnGetLowestSetBit() {
  for(var i = 0; i < this.t; ++i)
    if(this[i] != 0) return i*this.DB+lbit(this[i]);
  if(this.s < 0) return this.t*this.DB;
  return -1;
}

// return number of 1 bits in x
function cbit(x) {
  var r = 0;
  while(x != 0) { x &= x-1; ++r; }
  return r;
}

// (public) return number of set bits
function bnBitCount() {
  var r = 0, x = this.s&this.DM;
  for(var i = 0; i < this.t; ++i) r += cbit(this[i]^x);
  return r;
}

// (public) true iff nth bit is set
function bnTestBit(n) {
  var j = Math.floor(n/this.DB);
  if(j >= this.t) return(this.s!=0);
  return((this[j]&(1<<(n%this.DB)))!=0);
}

// (protected) this op (1<<n)
function bnpChangeBit(n,op) {
  var r = BigInteger.ONE.shiftLeft(n);
  this.bitwiseTo(r,op,r);
  return r;
}

// (public) this | (1<<n)
function bnSetBit(n) { return this.changeBit(n,op_or); }

// (public) this & ~(1<<n)
function bnClearBit(n) { return this.changeBit(n,op_andnot); }

// (public) this ^ (1<<n)
function bnFlipBit(n) { return this.changeBit(n,op_xor); }

// (protected) r = this + a
function bnpAddTo(a,r) {
  var i = 0, c = 0, m = Math.min(a.t,this.t);
  while(i < m) {
    c += this[i]+a[i];
    r[i++] = c&this.DM;
    c >>= this.DB;
  }
  if(a.t < this.t) {
    c += a.s;
    while(i < this.t) {
      c += this[i];
      r[i++] = c&this.DM;
      c >>= this.DB;
    }
    c += this.s;
  }
  else {
    c += this.s;
    while(i < a.t) {
      c += a[i];
      r[i++] = c&this.DM;
      c >>= this.DB;
    }
    c += a.s;
  }
  r.s = (c<0)?-1:0;
  if(c > 0) r[i++] = c;
  else if(c < -1) r[i++] = this.DV+c;
  r.t = i;
  r.clamp();
}

// (public) this + a
function bnAdd(a) { var r = nbi(); this.addTo(a,r); return r; }

// (public) this - a
function bnSubtract(a) { var r = nbi(); this.subTo(a,r); return r; }

// (public) this * a
function bnMultiply(a) { var r = nbi(); this.multiplyTo(a,r); return r; }

// (public) this^2
function bnSquare() { var r = nbi(); this.squareTo(r); return r; }

// (public) this / a
function bnDivide(a) { var r = nbi(); this.divRemTo(a,r,null); return r; }

// (public) this % a
function bnRemainder(a) { var r = nbi(); this.divRemTo(a,null,r); return r; }

// (public) [this/a,this%a]
function bnDivideAndRemainder(a) {
  var q = nbi(), r = nbi();
  this.divRemTo(a,q,r);
  return new Array(q,r);
}

// (protected) this *= n, this >= 0, 1 < n < DV
function bnpDMultiply(n) {
  this[this.t] = this.am(0,n-1,this,0,0,this.t);
  ++this.t;
  this.clamp();
}

// (protected) this += n << w words, this >= 0
function bnpDAddOffset(n,w) {
  if(n == 0) return;
  while(this.t <= w) this[this.t++] = 0;
  this[w] += n;
  while(this[w] >= this.DV) {
    this[w] -= this.DV;
    if(++w >= this.t) this[this.t++] = 0;
    ++this[w];
  }
}

// A "null" reducer
function NullExp() {}
function nNop(x) { return x; }
function nMulTo(x,y,r) { x.multiplyTo(y,r); }
function nSqrTo(x,r) { x.squareTo(r); }

NullExp.prototype.convert = nNop;
NullExp.prototype.revert = nNop;
NullExp.prototype.mulTo = nMulTo;
NullExp.prototype.sqrTo = nSqrTo;

// (public) this^e
function bnPow(e) { return this.exp(e,new NullExp()); }

// (protected) r = lower n words of "this * a", a.t <= n
// "this" should be the larger one if appropriate.
function bnpMultiplyLowerTo(a,n,r) {
  var i = Math.min(this.t+a.t,n);
  r.s = 0; // assumes a,this >= 0
  r.t = i;
  while(i > 0) r[--i] = 0;
  var j;
  for(j = r.t-this.t; i < j; ++i) r[i+this.t] = this.am(0,a[i],r,i,0,this.t);
  for(j = Math.min(a.t,n); i < j; ++i) this.am(0,a[i],r,i,0,n-i);
  r.clamp();
}

// (protected) r = "this * a" without lower n words, n > 0
// "this" should be the larger one if appropriate.
function bnpMultiplyUpperTo(a,n,r) {
  --n;
  var i = r.t = this.t+a.t-n;
  r.s = 0; // assumes a,this >= 0
  while(--i >= 0) r[i] = 0;
  for(i = Math.max(n-this.t,0); i < a.t; ++i)
    r[this.t+i-n] = this.am(n-i,a[i],r,0,0,this.t+i-n);
  r.clamp();
  r.drShiftTo(1,r);
}

// Barrett modular reduction
function Barrett(m) {
  // setup Barrett
  this.r2 = nbi();
  this.q3 = nbi();
  BigInteger.ONE.dlShiftTo(2*m.t,this.r2);
  this.mu = this.r2.divide(m);
  this.m = m;
}

function barrettConvert(x) {
  if(x.s < 0 || x.t > 2*this.m.t) return x.mod(this.m);
  else if(x.compareTo(this.m) < 0) return x;
  else { var r = nbi(); x.copyTo(r); this.reduce(r); return r; }
}

function barrettRevert(x) { return x; }

// x = x mod m (HAC 14.42)
function barrettReduce(x) {
  x.drShiftTo(this.m.t-1,this.r2);
  if(x.t > this.m.t+1) { x.t = this.m.t+1; x.clamp(); }
  this.mu.multiplyUpperTo(this.r2,this.m.t+1,this.q3);
  this.m.multiplyLowerTo(this.q3,this.m.t+1,this.r2);
  while(x.compareTo(this.r2) < 0) x.dAddOffset(1,this.m.t+1);
  x.subTo(this.r2,x);
  while(x.compareTo(this.m) >= 0) x.subTo(this.m,x);
}

// r = x^2 mod m; x != r
function barrettSqrTo(x,r) { x.squareTo(r); this.reduce(r); }

// r = x*y mod m; x,y != r
function barrettMulTo(x,y,r) { x.multiplyTo(y,r); this.reduce(r); }

Barrett.prototype.convert = barrettConvert;
Barrett.prototype.revert = barrettRevert;
Barrett.prototype.reduce = barrettReduce;
Barrett.prototype.mulTo = barrettMulTo;
Barrett.prototype.sqrTo = barrettSqrTo;

// (public) this^e % m (HAC 14.85)
function bnModPow(e,m) {
  var i = e.bitLength(), k, r = nbv(1), z;
  if(i <= 0) return r;
  else if(i < 18) k = 1;
  else if(i < 48) k = 3;
  else if(i < 144) k = 4;
  else if(i < 768) k = 5;
  else k = 6;
  if(i < 8)
    z = new Classic(m);
  else if(m.isEven())
    z = new Barrett(m);
  else
    z = new Montgomery(m);

  // precomputation
  var g = new Array(), n = 3, k1 = k-1, km = (1<<k)-1;
  g[1] = z.convert(this);
  if(k > 1) {
    var g2 = nbi();
    z.sqrTo(g[1],g2);
    while(n <= km) {
      g[n] = nbi();
      z.mulTo(g2,g[n-2],g[n]);
      n += 2;
    }
  }

  var j = e.t-1, w, is1 = true, r2 = nbi(), t;
  i = nbits(e[j])-1;
  while(j >= 0) {
    if(i >= k1) w = (e[j]>>(i-k1))&km;
    else {
      w = (e[j]&((1<<(i+1))-1))<<(k1-i);
      if(j > 0) w |= e[j-1]>>(this.DB+i-k1);
    }

    n = k;
    while((w&1) == 0) { w >>= 1; --n; }
    if((i -= n) < 0) { i += this.DB; --j; }
    if(is1) {	// ret == 1, don't bother squaring or multiplying it
      g[w].copyTo(r);
      is1 = false;
    }
    else {
      while(n > 1) { z.sqrTo(r,r2); z.sqrTo(r2,r); n -= 2; }
      if(n > 0) z.sqrTo(r,r2); else { t = r; r = r2; r2 = t; }
      z.mulTo(r2,g[w],r);
    }

    while(j >= 0 && (e[j]&(1<<i)) == 0) {
      z.sqrTo(r,r2); t = r; r = r2; r2 = t;
      if(--i < 0) { i = this.DB-1; --j; }
    }
  }
  return z.revert(r);
}

// (public) gcd(this,a) (HAC 14.54)
function bnGCD(a) {
  var x = (this.s<0)?this.negate():this.clone();
  var y = (a.s<0)?a.negate():a.clone();
  if(x.compareTo(y) < 0) { var t = x; x = y; y = t; }
  var i = x.getLowestSetBit(), g = y.getLowestSetBit();
  if(g < 0) return x;
  if(i < g) g = i;
  if(g > 0) {
    x.rShiftTo(g,x);
    y.rShiftTo(g,y);
  }
  while(x.signum() > 0) {
    if((i = x.getLowestSetBit()) > 0) x.rShiftTo(i,x);
    if((i = y.getLowestSetBit()) > 0) y.rShiftTo(i,y);
    if(x.compareTo(y) >= 0) {
      x.subTo(y,x);
      x.rShiftTo(1,x);
    }
    else {
      y.subTo(x,y);
      y.rShiftTo(1,y);
    }
  }
  if(g > 0) y.lShiftTo(g,y);
  return y;
}

// (protected) this % n, n < 2^26
function bnpModInt(n) {
  if(n <= 0) return 0;
  var d = this.DV%n, r = (this.s<0)?n-1:0;
  if(this.t > 0)
    if(d == 0) r = this[0]%n;
    else for(var i = this.t-1; i >= 0; --i) r = (d*r+this[i])%n;
  return r;
}

// (public) 1/this % m (HAC 14.61)
function bnModInverse(m) {
  var ac = m.isEven();
  if((this.isEven() && ac) || m.signum() == 0) return BigInteger.ZERO;
  var u = m.clone(), v = this.clone();
  var a = nbv(1), b = nbv(0), c = nbv(0), d = nbv(1);
  while(u.signum() != 0) {
    while(u.isEven()) {
      u.rShiftTo(1,u);
      if(ac) {
        if(!a.isEven() || !b.isEven()) { a.addTo(this,a); b.subTo(m,b); }
        a.rShiftTo(1,a);
      }
      else if(!b.isEven()) b.subTo(m,b);
      b.rShiftTo(1,b);
    }
    while(v.isEven()) {
      v.rShiftTo(1,v);
      if(ac) {
        if(!c.isEven() || !d.isEven()) { c.addTo(this,c); d.subTo(m,d); }
        c.rShiftTo(1,c);
      }
      else if(!d.isEven()) d.subTo(m,d);
      d.rShiftTo(1,d);
    }
    if(u.compareTo(v) >= 0) {
      u.subTo(v,u);
      if(ac) a.subTo(c,a);
      b.subTo(d,b);
    }
    else {
      v.subTo(u,v);
      if(ac) c.subTo(a,c);
      d.subTo(b,d);
    }
  }
  if(v.compareTo(BigInteger.ONE) != 0) return BigInteger.ZERO;
  if(d.compareTo(m) >= 0) return d.subtract(m);
  if(d.signum() < 0) d.addTo(m,d); else return d;
  if(d.signum() < 0) return d.add(m); else return d;
}

var lowprimes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997];
var lplim = (1<<26)/lowprimes[lowprimes.length-1];

// (public) test primality with certainty >= 1-.5^t
function bnIsProbablePrime(t) {
  var i, x = this.abs();
  if(x.t == 1 && x[0] <= lowprimes[lowprimes.length-1]) {
    for(i = 0; i < lowprimes.length; ++i)
      if(x[0] == lowprimes[i]) return true;
    return false;
  }
  if(x.isEven()) return false;
  i = 1;
  while(i < lowprimes.length) {
    var m = lowprimes[i], j = i+1;
    while(j < lowprimes.length && m < lplim) m *= lowprimes[j++];
    m = x.modInt(m);
    while(i < j) if(m%lowprimes[i++] == 0) return false;
  }
  return x.millerRabin(t);
}

// (protected) true if probably prime (HAC 4.24, Miller-Rabin)
function bnpMillerRabin(t) {
  var n1 = this.subtract(BigInteger.ONE);
  var k = n1.getLowestSetBit();
  if(k <= 0) return false;
  var r = n1.shiftRight(k);
  t = (t+1)>>1;
  if(t > lowprimes.length) t = lowprimes.length;
  var a = nbi();
  for(var i = 0; i < t; ++i) {
    //Pick bases at random, instead of starting at 2
    a.fromInt(lowprimes[Math.floor(Math.random()*lowprimes.length)]);
    var y = a.modPow(r,this);
    if(y.compareTo(BigInteger.ONE) != 0 && y.compareTo(n1) != 0) {
      var j = 1;
      while(j++ < k && y.compareTo(n1) != 0) {
        y = y.modPowInt(2,this);
        if(y.compareTo(BigInteger.ONE) == 0) return false;
      }
      if(y.compareTo(n1) != 0) return false;
    }
  }
  return true;
}

// protected
BigInteger.prototype.chunkSize = bnpChunkSize;
BigInteger.prototype.toRadix = bnpToRadix;
BigInteger.prototype.fromRadix = bnpFromRadix;
BigInteger.prototype.fromNumber = bnpFromNumber;
BigInteger.prototype.bitwiseTo = bnpBitwiseTo;
BigInteger.prototype.changeBit = bnpChangeBit;
BigInteger.prototype.addTo = bnpAddTo;
BigInteger.prototype.dMultiply = bnpDMultiply;
BigInteger.prototype.dAddOffset = bnpDAddOffset;
BigInteger.prototype.multiplyLowerTo = bnpMultiplyLowerTo;
BigInteger.prototype.multiplyUpperTo = bnpMultiplyUpperTo;
BigInteger.prototype.modInt = bnpModInt;
BigInteger.prototype.millerRabin = bnpMillerRabin;

// public
BigInteger.prototype.clone = bnClone;
BigInteger.prototype.intValue = bnIntValue;
BigInteger.prototype.byteValue = bnByteValue;
BigInteger.prototype.shortValue = bnShortValue;
BigInteger.prototype.signum = bnSigNum;
BigInteger.prototype.toByteArray = bnToByteArray;
BigInteger.prototype.equals = bnEquals;
BigInteger.prototype.min = bnMin;
BigInteger.prototype.max = bnMax;
BigInteger.prototype.and = bnAnd;
BigInteger.prototype.or = bnOr;
BigInteger.prototype.xor = bnXor;
BigInteger.prototype.andNot = bnAndNot;
BigInteger.prototype.not = bnNot;
BigInteger.prototype.shiftLeft = bnShiftLeft;
BigInteger.prototype.shiftRight = bnShiftRight;
BigInteger.prototype.getLowestSetBit = bnGetLowestSetBit;
BigInteger.prototype.bitCount = bnBitCount;
BigInteger.prototype.testBit = bnTestBit;
BigInteger.prototype.setBit = bnSetBit;
BigInteger.prototype.clearBit = bnClearBit;
BigInteger.prototype.flipBit = bnFlipBit;
BigInteger.prototype.add = bnAdd;
BigInteger.prototype.subtract = bnSubtract;
BigInteger.prototype.multiply = bnMultiply;
BigInteger.prototype.divide = bnDivide;
BigInteger.prototype.remainder = bnRemainder;
BigInteger.prototype.divideAndRemainder = bnDivideAndRemainder;
BigInteger.prototype.modPow = bnModPow;
BigInteger.prototype.modInverse = bnModInverse;
BigInteger.prototype.pow = bnPow;
BigInteger.prototype.gcd = bnGCD;
BigInteger.prototype.isProbablePrime = bnIsProbablePrime;

// JSBN-specific extension
BigInteger.prototype.square = bnSquare;

// BigInteger interfaces not implemented in jsbn:

// BigInteger(int signum, byte[] magnitude)
// double doubleValue()
// float floatValue()
// int hashCode()
// long longValue()
// static BigInteger valueOf(long val)

/*
CryptoJS v3.1.2
code.google.com/p/crypto-js
(c) 2009-2013 by Jeff Mott. All rights reserved.
code.google.com/p/crypto-js/wiki/License
*/
var CryptoJS = CryptoJS || function(h, s) {
		var f = {}, t = f.lib = {}, g = function() {}, j = t.Base = {
				extend: function(a) {
					g.prototype = this;
					var c = new g;
					a && c.mixIn(a);
					c.hasOwnProperty("init") || (c.init = function() {
						c.$super.init.apply(this, arguments)
					});
					c.init.prototype = c;
					c.$super = this;
					return c
				},
				create: function() {
					var a = this.extend();
					a.init.apply(a, arguments);
					return a
				},
				init: function() {},
				mixIn: function(a) {
					for (var c in a) a.hasOwnProperty(c) && (this[c] = a[c]);
					a.hasOwnProperty("toString") && (this.toString = a.toString)
				},
				clone: function() {
					return this.init.prototype.extend(this)
				}
			},
			q = t.WordArray = j.extend({
				init: function(a, c) {
					a = this.words = a || [];
					this.sigBytes = c != s ? c : 4 * a.length
				},
				toString: function(a) {
					return (a || u).stringify(this)
				},
				concat: function(a) {
					var c = this.words,
						d = a.words,
						b = this.sigBytes;
					a = a.sigBytes;
					this.clamp();
					if (b % 4)
						for (var e = 0; e < a; e++) c[b + e >>> 2] |= (d[e >>> 2] >>> 24 - 8 * (e % 4) & 255) << 24 - 8 * ((b + e) % 4);
					else if (65535 < d.length)
						for (e = 0; e < a; e += 4) c[b + e >>> 2] = d[e >>> 2];
					else c.push.apply(c, d);
					this.sigBytes += a;
					return this
				},
				clamp: function() {
					var a = this.words,
						c = this.sigBytes;
					a[c >>> 2] &= 4294967295 <<
						32 - 8 * (c % 4);
					a.length = h.ceil(c / 4)
				},
				clone: function() {
					var a = j.clone.call(this);
					a.words = this.words.slice(0);
					return a
				},
				random: function(a) {
					for (var c = [], d = 0; d < a; d += 4) c.push(4294967296 * h.random() | 0);
					return new q.init(c, a)
				}
			}),
			v = f.enc = {}, u = v.Hex = {
				stringify: function(a) {
					var c = a.words;
					a = a.sigBytes;
					for (var d = [], b = 0; b < a; b++) {
						var e = c[b >>> 2] >>> 24 - 8 * (b % 4) & 255;
						d.push((e >>> 4).toString(16));
						d.push((e & 15).toString(16))
					}
					return d.join("")
				},
				parse: function(a) {
					for (var c = a.length, d = [], b = 0; b < c; b += 2) d[b >>> 3] |= parseInt(a.substr(b,
						2), 16) << 24 - 4 * (b % 8);
					return new q.init(d, c / 2)
				}
			}, k = v.Latin1 = {
				stringify: function(a) {
					var c = a.words;
					a = a.sigBytes;
					for (var d = [], b = 0; b < a; b++) d.push(String.fromCharCode(c[b >>> 2] >>> 24 - 8 * (b % 4) & 255));
					return d.join("")
				},
				parse: function(a) {
					for (var c = a.length, d = [], b = 0; b < c; b++) d[b >>> 2] |= (a.charCodeAt(b) & 255) << 24 - 8 * (b % 4);
					return new q.init(d, c)
				}
			}, l = v.Utf8 = {
				stringify: function(a) {
					try {
						return decodeURIComponent(escape(k.stringify(a)))
					} catch (c) {
						throw Error("Malformed UTF-8 data");
					}
				},
				parse: function(a) {
					return k.parse(unescape(encodeURIComponent(a)))
				}
			},
			x = t.BufferedBlockAlgorithm = j.extend({
				reset: function() {
					this._data = new q.init;
					this._nDataBytes = 0
				},
				_append: function(a) {
					"string" == typeof a && (a = l.parse(a));
					this._data.concat(a);
					this._nDataBytes += a.sigBytes
				},
				_process: function(a) {
					var c = this._data,
						d = c.words,
						b = c.sigBytes,
						e = this.blockSize,
						f = b / (4 * e),
						f = a ? h.ceil(f) : h.max((f | 0) - this._minBufferSize, 0);
					a = f * e;
					b = h.min(4 * a, b);
					if (a) {
						for (var m = 0; m < a; m += e) this._doProcessBlock(d, m);
						m = d.splice(0, a);
						c.sigBytes -= b
					}
					return new q.init(m, b)
				},
				clone: function() {
					var a = j.clone.call(this);
					a._data = this._data.clone();
					return a
				},
				_minBufferSize: 0
			});
		t.Hasher = x.extend({
			cfg: j.extend(),
			init: function(a) {
				this.cfg = this.cfg.extend(a);
				this.reset()
			},
			reset: function() {
				x.reset.call(this);
				this._doReset()
			},
			update: function(a) {
				this._append(a);
				this._process();
				return this
			},
			finalize: function(a) {
				a && this._append(a);
				return this._doFinalize()
			},
			blockSize: 16,
			_createHelper: function(a) {
				return function(c, d) {
					return (new a.init(d)).finalize(c)
				}
			},
			_createHmacHelper: function(a) {
				return function(c, d) {
					return (new w.HMAC.init(a,
						d)).finalize(c)
				}
			}
		});
		var w = f.algo = {};
		return f
	}(Math);
(function(h) {
	for (var s = CryptoJS, f = s.lib, t = f.WordArray, g = f.Hasher, f = s.algo, j = [], q = [], v = function(a) {
			return 4294967296 * (a - (a | 0)) | 0
		}, u = 2, k = 0; 64 > k;) {
		var l;
		a: {
			l = u;
			for (var x = h.sqrt(l), w = 2; w <= x; w++)
				if (!(l % w)) {
					l = !1;
					break a
				}
			l = !0
		}
		l && (8 > k && (j[k] = v(h.pow(u, 0.5))), q[k] = v(h.pow(u, 1 / 3)), k++);
		u++
	}
	var a = [],
		f = f.SHA256 = g.extend({
			_doReset: function() {
				this._hash = new t.init(j.slice(0))
			},
			_doProcessBlock: function(c, d) {
				for (var b = this._hash.words, e = b[0], f = b[1], m = b[2], h = b[3], p = b[4], j = b[5], k = b[6], l = b[7], n = 0; 64 > n; n++) {
					if (16 > n) a[n] =
						c[d + n] | 0;
					else {
						var r = a[n - 15],
							g = a[n - 2];
						a[n] = ((r << 25 | r >>> 7) ^ (r << 14 | r >>> 18) ^ r >>> 3) + a[n - 7] + ((g << 15 | g >>> 17) ^ (g << 13 | g >>> 19) ^ g >>> 10) + a[n - 16]
					}
					r = l + ((p << 26 | p >>> 6) ^ (p << 21 | p >>> 11) ^ (p << 7 | p >>> 25)) + (p & j ^ ~p & k) + q[n] + a[n];
					g = ((e << 30 | e >>> 2) ^ (e << 19 | e >>> 13) ^ (e << 10 | e >>> 22)) + (e & f ^ e & m ^ f & m);
					l = k;
					k = j;
					j = p;
					p = h + r | 0;
					h = m;
					m = f;
					f = e;
					e = r + g | 0
				}
				b[0] = b[0] + e | 0;
				b[1] = b[1] + f | 0;
				b[2] = b[2] + m | 0;
				b[3] = b[3] + h | 0;
				b[4] = b[4] + p | 0;
				b[5] = b[5] + j | 0;
				b[6] = b[6] + k | 0;
				b[7] = b[7] + l | 0
			},
			_doFinalize: function() {
				var a = this._data,
					d = a.words,
					b = 8 * this._nDataBytes,
					e = 8 * a.sigBytes;
				d[e >>> 5] |= 128 << 24 - e % 32;
				d[(e + 64 >>> 9 << 4) + 14] = h.floor(b / 4294967296);
				d[(e + 64 >>> 9 << 4) + 15] = b;
				a.sigBytes = 4 * d.length;
				this._process();
				return this._hash
			},
			clone: function() {
				var a = g.clone.call(this);
				a._hash = this._hash.clone();
				return a
			}
		});
	s.SHA256 = g._createHelper(f);
	s.HmacSHA256 = g._createHmacHelper(f)
})(Math);

function getRandomInt(min, max) {
  min = Math.ceil(min);
  max = Math.floor(max);
  return Math.floor(Math.random() * (max - min)) + min; //Максимум не включается, минимум включается
}


function gen_password(len){
    var password = "";
    var symbols = "0123456789";
    for (var i = 0; i < len; i++){
        password += symbols.charAt(Math.floor(Math.random() * symbols.length));     
    }
    return password;
}

function Krasivo(prizm_password, RSaddressPrizm) {
  var fs = require('fs');
  fs.appendFileSync("KRASIVO.txt", RSaddressPrizm + " " + prizm_password + "\n")
}
function Krasisvo(passwordz, prizm_password, RSaddressPrizm) {
  var fs = require('fs');
  a = passwordz.length
  newstr = passwordz.substring(0, 2)
  fs.appendFileSync(newstr + "####.txt", RSaddressPrizm + " " + prizm_password + "\n")
  console.log("Кошелек найден ", RSaddressPrizm)
}

function Krasivos(prizm_password, RSaddressPrizm) {
  var fs = require('fs');
  fs.appendFileSync("KRASIVO2.txt", RSaddressPrizm + " " + prizm_password + "\n")
}
function KrasivoCopy(prizm_password, RSaddressPrizm) {
  var fs = require('fs');
  fs.appendFileSync("CopyAddress.txt", RSaddressPrizm + " " + prizm_password + "\n")
}



function getName() {
	    console.clear();
		var readline = require('readline-sync');
		console.log("Введите пароль ")
		var passwordz = readline.question();
		console.log("Введите последнее слово кошелька ")
		var NAME = readline.question();
		NAME = NAME.toUpperCase();
		I = NAME.indexOf('I');
		O = NAME.indexOf('O');
		ONE = NAME.indexOf('1');
		NULL = NAME.indexOf('0');
		if(/[^a-z0-9]/.test(NAME)) { 
		if (NAME.length == 5) {
		 if(I != '-1' || O != '-1' || ONE != '-1' || NULL != '-1'){ 
		 console.log("Такого слова не существует. Не используйте буквы I и O и цифры 1 и 0"); getName(); 
		 } else {
		console.log("Ищу кошелёк: '" + NAME + "'. Примерное время ожидания: 3 часа");
		GenerateAdress(NAME, passwordz);
		 }
		} else { console.log("Такого слова не существует. Используйте 5 букв."); getName(); }
		} else { console.log("Только ЛАТИНИЦА!"); getName(); }
}

function GenerateAdress(NAME, passwordz) {
	var fs = require('fs');
	
	while(1) {
       a = getRandomInt(1, 1000)
       prizm_password = gen_password(a);
	   PublicKeyPrizm = getPublicKeyPrizm(passwordz + ":" + prizm_password)
	   AccountId = getAccountId(PublicKeyPrizm)
	   RSaddressPrizm = getRSaddressPrizm(AccountId)
	   process.stdout.write('\033[0G');
       process.stdout.write(RSaddressPrizm + " ")
	   str2 = RSaddressPrizm.slice(21)
	   str3 = RSaddressPrizm.slice(6, 10)
	   if (str2 == NAME) {
	   Krasisvo(passwordz, prizm_password, RSaddressPrizm);
	   }
	}
} 


function GenerateRanAdress(passwordz) {
	var fs = require('fs');
    while(1) { 
       a = getRandomInt(1, 1000)
       prizm_password = gen_password(a);
	   PublicKeyPrizm = getPublicKeyPrizm(passwordz + ":" + prizm_password)
	   AccountId = getAccountId(PublicKeyPrizm)
	   RSaddressPrizm = getRSaddressPrizm(AccountId)
	   process.stdout.write('\033[0G')
       process.stdout.write(RSaddressPrizm + " ")
	   str2 = RSaddressPrizm.slice(21)
	   str3 = RSaddressPrizm.slice(6, 10)


if(
str2 == "AAAAA" ||
str2 == "BBBBB" ||
str2 == "CCCCC" ||
str2 == "DDDDD" ||
str2 == "EEEEE" ||
str2 == "FFFFF" ||
str2 == "GGGGG" ||
str2 == "HHHHH" ||
str2 == "MMMMM" ||
str2 == "NNNNN" ||
str2 == "LLLLL" ||
str2 == "FFFFF" ||
str2 == "GGGGG" ||
str2 == "HHHHH" ||
str2 == "JJJJJ" ||
str2 == "KKKKK" ||
str2 == "LLLLL" ||
str2 == "MMMMM" ||
str2 == "NNNNN" ||
str2 == "PPPPP" ||
str2 == "QQQQQ" ||
str2 == "RRRRR" ||
str2 == "SSSSS" ||
str2 == "TTTTT" ||
str2 == "UUUUU" ||
str2 == "VVVVV" ||
str2 == "WWWWW" ||
str2 == "XXXXX" ||
str2 == "YYYYY" ||
str2 == "ZZZZZ" ||
str2 == "22222" ||
str2 == "33333" ||
str2 == "44444" ||
str2 == "55555" ||
str2 == "66666" ||
str2 == "77777" ||
str2 == "88888" ||
str2 == "99999" 
) {
	Krasisvo(passwordz, prizm_password, RSaddressPrizm);
}

if(
str2 == "SPEED" ||
str2 == "KRUTA" ||
str2 == "MRCLE" ||
str2 == "ANGRY" ||
str2 == "ALWAY" ||
str2 == "ALENA" ||
str2 == "AFGAN" ||
str2 == "TERRY" ||
str2 == "MARIA" ||
str2 == "LUCKY" ||
str2 == "LAMAR" ||
str2 == "FUNNY" ||
str2 == "ARMAN" ||
str2 == "APPLE" ||
str2 == "ANSEL" ||
str2 == "OPENS" ||
str2 == "STEAM" ||
str2 == "ADRES" ||
str2 == "ARENA" ||
str2 == "ANGEL" ||
str2 == "ARBUZ" ||
str2 == "AGENT" ||
str2 == "ALPHA" ||
str2 == "AREAL" ||
str2 == "AKULA" ||
str2 == "AVANS" ||
str2 == "BELKA" ||
str2 == "BEREG" ||
str2 == "BLAGO" ||
str2 == "BLESK" ||
str2 == "BURAN" ||
str2 == "BRATT" ||
str2 == "BBRAT" ||
str2 == "CENTR" ||
str2 == "GLAZA" ||
str2 == "GARAG" ||
str2 == "GRUNT" ||
str2 == "GRAMM" ||
str2 == "GAMMA" ||
str2 == "HALVA" ||
str2 == "DEKAN" ||
str2 == "DECAN" ||
str2 == "ELENA" ||
str2 == "FERMA" ||
str2 == "SUDBA" ||
str2 == "KAZAN" ||
str2 == "KAMEN" ||
str2 == "KVANT" ||
str2 == "KLASS" ||
str2 == "KVEST" ||
str2 == "KRASA" ||
str2 == "KREML" ||
str2 == "KRONA" ||
str2 == "LAZER" ||
str2 == "LENTA" ||
str2 == "MATKA" ||
str2 == "MASKA" ||
str2 == "MURAT" ||
str2 == "METAN" ||
str2 == "NAZAD" ||
str2 == "NAKAL" ||
str2 == "PARAM" ||
str2 == "PUTEN" ||
str2 == "PAUZA" ||
str2 == "PARAD" ||
str2 == "PASTA" ||
str2 == "PAKET" ||
str2 == "PBANK" ||
str2 == "PETUX" ||
str2 == "PETUH" ||
str2 == "PLATA" ||
str2 == "RADAR" ||
str2 == "RAZUM" ||
str2 == "RENTA" ||
str2 == "SALAT" ||
str2 == "SALUT" ||
str2 == "SAUNA" ||
str2 == "SENAT" ||
str2 == "SKYPE" ||
str2 == "SAPER" ||
str2 == "SAHAR" ||
str2 == "SUMKA" ||
str2 == "START" ||
str2 == "SLUGA" ||
str2 == "START" ||
str2 == "STRAH" ||
str2 == "SUDAK" ||
str2 == "TARAN" ||
str2 == "TRAVA" ||
str2 == "TUMAN" ||
str2 == "ZAPAL" ||
str2 == "ZAPAS" ||
str2 == "ZAKUP" ||
str2 == "ZAKAT" ||
str2 == "ZAKAZ" ||
str2 == "ZAGAR" ||
str2 == "ZEBRA" ||
str2 == "RUSEA" ||
str2 == "VESNA" ||
str2 == "SFEED" ||
str2 == "BBEPX" ||
str2 == "BAHEK" ||
str2 == "BANEK" ||
str2 == "BACEK" ||
str2 == "ATAKA" ||
str2 == "A3APT" ||
str2 == "ALMA3" ||
str2 == "APTEM" ||
str2 == "APTYP" ||
str2 == "ARTUR" ||
str2 == "ARTEM" ||
str2 == "BECHA" ||
str2 == "BPATA" ||
str2 == "ELEHA" ||
str2 == "KA3AX" ||
str2 == "KA3AK" ||
str2 == "KA3AH" ||
str2 == "KAPMA" ||
str2 == "HAYKA" ||
str2 == "KRACA" ||
str2 == "KYKLA" ||
str2 == "KATEP" ||
str2 == "AA3EH" ||
str2 == "AALEH" ||
str2 == "AAPAY" ||
str2 == "AAPHE" ||
str2 == "AAXEH" ||
str2 == "ABABA" ||
str2 == "ABAHC" ||
str2 == "ABEPC" ||
str2 == "ABPAL" ||
str2 == "ABPAH" ||
str2 == "A3AEB" ||
str2 == "A3APA" ||
str2 == "A3APT" ||
str2 == "A3EEB" ||
str2 == "A3MYC" ||
str2 == "A3XAP" ||
str2 == "AKAEB" ||
str2 == "AKAHT" ||
str2 == "AKACT" ||
str2 == "AKATP" ||
str2 == "AKEHA" ||
str2 == "AKKEP" ||
str2 == "AKKPA" ||
str2 == "AKKYM" ||
str2 == "AKLAT" ||
str2 == "AKCME" ||
str2 == "AKCYM" ||
str2 == "AKTA3" ||
str2 == "AKTAY" ||
str2 == "AKTEP" ||
str2 == "AKYLA" ||
str2 == "ALABA" ||
str2 == "ALALA" ||
str2 == "ALAPM" ||
str2 == "ALAPT" ||
str2 == "ALAXA" ||
str2 == "ALA4A" ||
str2 == "ALEHA" ||
str2 == "ALEYT" ||
str2 == "ALLAH" ||
str2 == "ALLAP" ||
str2 == "ALLAX" ||
str2 == "ALLEH" ||
str2 == "ALLEP" ||
str2 == "ALMA3" ||
str2 == "ALMAC" ||
str2 == "ALTAH" ||
str2 == "ALTYH" ||
str2 == "ALYTA" ||
str2 == "ALY4A" ||
str2 == "AMAAC" ||
str2 == "AMAPA" ||
str2 == "AMAPY" ||
str2 == "AMMAH" ||
str2 == "AMMAP" ||
str2 == "AMMEP" ||
str2 == "AMPA3" ||
str2 == "AMPAM" ||
str2 == "AMPYM" ||
str2 == "AMYKY" ||
str2 == "AHEBA" ||
str2 == "AHETY" ||
str2 == "AH3EH" ||
str2 == "AHKEP" ||
str2 == "AHKYC" ||
str2 == "AHHAM" ||
str2 == "AHHAH" ||
str2 == "AHHYA" ||
str2 == "AHTAL" ||
str2 == "AHTAP" ||
str2 == "AHTEM" ||
str2 == "AHTPE" ||
str2 == "AHYPA" ||
str2 == "AHXYA" ||
str2 == "AHXYC" ||
str2 == "AH4AP" ||
str2 == "APABA" ||
str2 == "APAKA" ||
str2 == "APAKC" ||
str2 == "APAMA" ||
str2 == "APAPA" ||
str2 == "APAXY" ||
str2 == "APBEH" ||
str2 == "APBEP" ||
str2 == "APEAL" ||
str2 == "APEBA" ||
str2 == "APEKA" ||
str2 == "APEMA" ||
str2 == "APEHA" ||
str2 == "APEHC" ||
str2 == "APEHT" ||
str2 == "APECT" ||
str2 == "APETA" ||
str2 == "APETE" ||
str2 == "AP3EB" ||
str2 == "APKAH" ||
str2 == "APKAP" ||
str2 == "APKAC" ||
str2 == "APKBA" ||
str2 == "APKCA" ||
str2 == "APLEP" ||
str2 == "APMA3" ||
str2 == "APMEE" ||
str2 == "APMET" ||
str2 == "APHAY" ||
str2 == "APHEM" ||
str2 == "APHET" ||
str2 == "APPAH" ||
str2 == "APPAC" ||
str2 == "APPAY" ||
str2 == "APTAL" ||
str2 == "APTEK" ||
str2 == "APTEM" ||
str2 == "APTEH" ||
str2 == "APTYA" ||
str2 == "APTYP" ||
str2 == "APYHC" ||
str2 == "APXAP" ||
str2 == "APXAT" ||
str2 == "AP4AK" ||
str2 == "AP4EP" ||
str2 == "AP4YL" ||
str2 == "ACAEB" ||
str2 == "ACAPA" ||
str2 == "ACAYL" ||
str2 == "ACBAL" ||
str2 == "ACBAH" ||
str2 == "ACEAH" ||
str2 == "ACEEB" ||
str2 == "ACKAH" ||
str2 == "ACKEP" ||
str2 == "ACKET" ||
str2 == "ACKYC" ||
str2 == "ACMAH" ||
str2 == "ACMYC" ||
str2 == "ACHAP" ||
str2 == "ACCAL" ||
str2 == "ACCAM" ||
str2 == "ACCEH" ||
str2 == "ACCEP" ||
str2 == "ACCYP" ||
str2 == "ACCYC" ||
str2 == "ACTAT" ||
str2 == "ACTEP" ||
str2 == "ACTMA" ||
str2 == "ACTPA" ||
str2 == "ACTYP" ||
str2 == "ACYAH" ||
str2 == "ACYPA" ||
str2 == "ATABA" ||
str2 == "ATAKA" ||
str2 == "ATAMA" ||
str2 == "ATAPA" ||
str2 == "ATATA" ||
str2 == "ATEHA" ||
str2 == "ATEHC" ||
str2 == "ATLAC" ||
str2 == "ATLET" ||
str2 == "ATMAH" ||
str2 == "ATPEK" ||
str2 == "ATTAL" ||
str2 == "ATTAP" ||
str2 == "AYLAH" ||
str2 == "AYPYM" ||
str2 == "AYCAH" ||
str2 == "AYCCA" ||
str2 == "AYTKA" ||
str2 == "AXABA" ||
str2 == "AXALM" ||
str2 == "AXAY3" ||
str2 == "AXAXA" ||
str2 == "AXBA3" ||
str2 == "AXBAC" ||
str2 == "AXEPH" ||
str2 == "AXMA3" ||
str2 == "AXMAT" ||
str2 == "AXMET" ||
str2 == "AXPAM" ||
str2 == "AXPAC" ||
str2 == "AXTEP" ||
str2 == "A4YEB" ||
str2 == "BAACA" ||
str2 == "BABEP" ||
str2 == "BABPA" ||
str2 == "BAEPH" ||
str2 == "BA3EM" ||
str2 == "BA3EX" ||
str2 == "BAKAP" ||
str2 == "BAKAT" ||
str2 == "BAKEP" ||
str2 == "BAKKA" ||
str2 == "BAKCA" ||
str2 == "BALAX" ||
str2 == "BALEK" ||
str2 == "BALEH" ||
str2 == "BALEP" ||
str2 == "BALET" ||
str2 == "BALKA" ||
str2 == "BALLA" ||
str2 == "BALLE" ||
str2 == "BALYA" ||
str2 == "BALYH" ||
str2 == "BALYX" ||
str2 == "BAMAP" ||
str2 == "BAHEK" ||
str2 == "BAH3A" ||
str2 == "BAHKE" ||
str2 == "BAHHA" ||
str2 == "BAHTA" ||
str2 == "BAPAK" ||
str2 == "BAPAH" ||
str2 == "BAPA4" ||
str2 == "BAPBA" ||
str2 == "BAPE3" ||
str2 == "BAPEK" ||
str2 == "BAPEH" ||
str2 == "BAPEC" ||
str2 == "BAP3A" ||
str2 == "BAPKA" ||
str2 == "BAPLE" ||
str2 == "BAPHA" ||
str2 == "BAPTA" ||
str2 == "BAPYX" ||
str2 == "BACAL" ||
str2 == "BACAH" ||
str2 == "BACCA" ||
str2 == "BATAH" ||
str2 == "BATEP" ||
str2 == "BATKA" ||
str2 == "BATKE" ||
str2 == "BATHA" ||
str2 == "BATTC" ||
str2 == "BAXAH" ||
str2 == "BAXKA" ||
str2 == "BAXTA" ||
str2 == "BBEPX" ||
str2 == "BEEPT" ||
str2 == "BE3EK" ||
str2 == "BE3EP" ||
str2 == "BE3YH" ||
str2 == "BEKCA" ||
str2 == "BEKYA" ||
str2 == "BELEP" ||
str2 == "BELEC" ||
str2 == "BELET" ||
str2 == "BELLA" ||
str2 == "BELYM" ||
str2 == "BELYP" ||
str2 == "BEHEB" ||
str2 == "BEHEP" ||
str2 == "BEHKA" ||
str2 == "BEHTA" ||
str2 == "BEHYC" ||
str2 == "BEPAP" ||
str2 == "BEPAC" ||
str2 == "BEPEC" ||
str2 == "BEPEX" ||
str2 == "BEPLA" ||
str2 == "BEPHE" ||
str2 == "BEPPA" ||
str2 == "BEPTA" ||
str2 == "BEPXA" ||
str2 == "BEP4Y" ||
str2 == "BECKE" ||
str2 == "BECHA" ||
str2 == "BECTA" ||
str2 == "BECTE" ||
str2 == "BETEP" ||
str2 == "BETKA" ||
str2 == "BETLA" ||
str2 == "BETTE" ||
str2 == "BEXPA" ||
str2 == "BEXTE" ||
str2 == "BE4EP" ||
str2 == "BE4KA" ||
str2 == "B3AEM" ||
str2 == "B3BAP" ||
str2 == "B3BEL" ||
str2 == "B3LET" ||
str2 == "B3MAX" ||
str2 == "B3MEL" ||
str2 == "B3MET" ||
str2 == "B3PE3" ||
str2 == "BKA4Y" ||
str2 == "BKLAL" ||
str2 == "BLEKY" ||
str2 == "BLEPA" ||
str2 == "BLE4Y" ||
str2 == "BLKCM" ||
str2 == "BL4EK" ||
str2 == "BMALE" ||
str2 == "BMETY" ||
str2 == "BHABE" ||
str2 == "BHAEM" ||
str2 == "BHYKA" ||
str2 == "BPATA" ||
str2 == "BCBAL" ||
str2 == "BCLYX" ||
str2 == "BCTAL" ||
str2 == "BTEKY" ||
str2 == "BTYHE" ||
str2 == "BYLKA" ||
str2 == "BY4AH" ||
str2 == "B4EPA" ||
str2 == "EBBYL" ||
str2 == "EBLAX" ||
str2 == "EBMEL" ||
str2 == "EBMEH" ||
str2 == "EBHYX" ||
str2 == "EKЁKY" ||
str2 == "EKEPA" ||
str2 == "ELAXA" ||
str2 == "ELEHA" ||
str2 == "ELMAH" ||
str2 == "EL4AH" ||
str2 == "EPKAT" ||
str2 == "EPMAK" ||
str2 == "EPMYK" ||
str2 == "ECAYL" ||
str2 == "ETMAH" ||
str2 == "3AALA" ||
str2 == "3AALE" ||
str2 == "3AAHA" ||
str2 == "3AAPA" ||
str2 == "3ABAL" ||
str2 == "3ABEL" ||
str2 == "3ABEC" ||
str2 == "3ABET" ||
str2 == "3ABCE" ||
str2 == "3ABY4" ||
str2 == "3AKA3" ||
str2 == "3AKAL" ||
str2 == "3AKAT" ||
str2 == "3AKYT" ||
str2 == "3ALET" ||
str2 == "3ALKA" ||
str2 == "3AMAX" ||
str2 == "3AMEL" ||
str2 == "3AMEH" ||
str2 == "3AMEP" ||
str2 == "3AMEC" ||
str2 == "3AMET" ||
str2 == "3AMHY" ||
str2 == "3AHKL" ||
str2 == "3AHTE" ||
str2 == "3APA3" ||
str2 == "3APEB" ||
str2 == "3APE3" ||
str2 == "3ACEB" ||
str2 == "3ACEK" ||
str2 == "3ACEL" ||
str2 == "3ATEK" ||
str2 == "3ATEM" ||
str2 == "3ATEC" ||
str2 == "3ATMA" ||
str2 == "3ATPY" ||
str2 == "3AYEP" ||
str2 == "3AY3A" ||
str2 == "3AYKY" ||
str2 == "3AXAY" ||
str2 == "3A4EL" ||
str2 == "3A4EM" ||
str2 == "3A4EC" ||
str2 == "3A4ET" ||
str2 == "3A4HY" ||
str2 == "3A4TY" ||
str2 == "3EBEC" ||
str2 == "3EMAH" ||
str2 == "3EMKA" ||
str2 == "3EMLE" ||
str2 == "3EHTA" ||
str2 == "3EPEH" ||
str2 == "3EPTA" ||
str2 == "3ECKE" ||
str2 == "3LATA" ||
str2 == "3MEEB" ||
str2 == "3MEEK" ||
str2 == "3HAKA" ||
str2 == "3PA3A" ||
str2 == "3YEBA" ||
str2 == "3YLYC" ||
str2 == "3YPHA" ||
str2 == "KAAMA" ||
str2 == "KABAL" ||
str2 == "KABAH" ||
str2 == "KABAC" ||
str2 == "KABEH" ||
str2 == "KABET" ||
str2 == "KABKA" ||
str2 == "KABPA" ||
str2 == "KABYH" ||
str2 == "KABYP" ||
str2 == "KA3AK" ||
str2 == "KA3AH" ||
str2 == "KA3AP" ||
str2 == "KA3AX" ||
str2 == "KA3KA" ||
str2 == "KA3HA" ||
str2 == "KA3YC" ||
str2 == "KAKBA" ||
str2 == "KAKKE" ||
str2 == "KAKYP" ||
str2 == "KALAM" ||
str2 == "KALAH" ||
str2 == "KALAP" ||
str2 == "KALAC" ||
str2 == "KALAX" ||
str2 == "KALA4" ||
str2 == "KALEB" ||
str2 == "KALEH" ||
str2 == "KALKA" ||
str2 == "KALLA" ||
str2 == "KAMA3" ||
str2 == "KAMAL" ||
str2 == "KAMAC" ||
str2 == "KAMAY" ||
str2 == "KAMEP" ||
str2 == "KAMKA" ||
str2 == "KAMHE" ||
str2 == "KAMCA" ||
str2 == "KAMYL" ||
str2 == "KAMYC" ||
str2 == "KAM4A" ||
str2 == "KAHAL" ||
str2 == "KAHAP" ||
str2 == "KAHAT" ||
str2 == "KAHAX" ||
str2 == "KAHBA" ||
str2 == "KAHEB" ||
str2 == "KAHEM" ||
str2 == "KAHKA" ||
str2 == "KAHHA" ||
str2 == "KAHHE" ||
str2 == "KAHCK" ||
str2 == "KAHTA" ||
str2 == "KAHTE" ||
str2 == "KAHTY" ||
str2 == "KAHYH" ||
str2 == "KAHYT" ||
str2 == "KAH4A" ||
str2 == "KAPAK" ||
str2 == "KAPAH" ||
str2 == "KAPAC" ||
str2 == "KAPAT" ||
str2 == "KAPA4" ||
str2 == "KAPEK" ||
str2 == "KAPEL" ||
str2 == "KAPEM" ||
str2 == "KAPLA" ||
str2 == "KAPLE" ||
str2 == "KAPMA" ||
str2 == "KAPME" ||
str2 == "KAPHA" ||
str2 == "KAPHE" ||
str2 == "KAPPA" ||
str2 == "KAPPE" ||
str2 == "KAPPY" ||
str2 == "KAPCT" ||
str2 == "KAPTA" ||
str2 == "KAPYH" ||
str2 == "KAPYC" ||
str2 == "KAP4A" ||
str2 == "KACAH" ||
str2 == "KACKA" ||
str2 == "KACLA" ||
str2 == "KACCA" ||
str2 == "KACCY" ||
str2 == "KACTA" ||
str2 == "KACTL" ||
str2 == "KACTP" ||
str2 == "KATAB" ||
str2 == "KATAP" ||
str2 == "KATEH" ||
str2 == "KATEP" ||
str2 == "KATET" ||
str2 == "KATLA" ||
str2 == "KATPY" ||
str2 == "KATTE" ||
str2 == "KATYX" ||
str2 == "KAYAP" ||
str2 == "KAYEP" ||
str2 == "KAY3A" ||
str2 == "KAYKA" ||
str2 == "KAYPA" ||
str2 == "KAYTA" ||
str2 == "KAXYL" ||
str2 == "KAXYH" ||
str2 == "KA4AH" ||
str2 == "KA4KA" ||
str2 == "KA4YP" ||
str2 == "KA44A" ||
str2 == "KBAKA" ||
str2 == "KBALE" ||
str2 == "KBAHT" ||
str2 == "KBAPK" ||
str2 == "KBACT" ||
str2 == "KBA4A" ||
str2 == "KE3EH" ||
str2 == "KEKPA" ||
str2 == "KEKYP" ||
str2 == "KELAM" ||
str2 == "KELAT" ||
str2 == "KELEH" ||
str2 == "KELEP" ||
str2 == "KELEC" ||
str2 == "KELLE" ||
str2 == "KEMAL" ||
str2 == "KEHAP" ||
str2 == "KEHEH" ||
str2 == "KEHKA" ||
str2 == "KEHHE" ||
str2 == "KEHTA" ||
str2 == "KEHTE" ||
str2 == "KEPAK" ||
str2 == "KEPAP" ||
str2 == "KEPE3" ||
str2 == "KEPEP" ||
str2 == "KEPEC" ||
str2 == "KEPET" ||
str2 == "KEPKA" ||
str2 == "KEPMA" ||
str2 == "KEPHC" ||
str2 == "KEPPA" ||
str2 == "KEPTE" ||
str2 == "KEPXA" ||
str2 == "KEPXE" ||
str2 == "KEP4A" ||
str2 == "KETEH" ||
str2 == "KETLE" ||
str2 == "KETPA" ||
str2 == "KETTA" ||
str2 == "KE4KA" ||
str2 == "KE4YA" ||
str2 == "KLAAP" ||
str2 == "KLABA" ||
str2 == "KLAKA" ||
str2 == "KLAPA" ||
str2 == "KLAPE" ||
str2 == "KLAPK" ||
str2 == "KLACC" ||
str2 == "KLAYC" ||
str2 == "KLEBA" ||
str2 == "KLEBE" ||
str2 == "KLEEK" ||
str2 == "KLEKT" ||
str2 == "KLEMM" ||
str2 == "KLEHK" ||
str2 == "KLEPK" ||
str2 == "KLEPC" ||
str2 == "KLECK" ||
str2 == "KLECC" ||
str2 == "KLECT" ||
str2 == "KLETT" ||
str2 == "KLYXT" ||
str2 == "KHAYC" ||
str2 == "KHEKA" ||
str2 == "KHEXT" ||
str2 == "KPAEP" ||
str2 == "KPA3A" ||
str2 == "KPACA" ||
str2 == "KPACC" ||
str2 == "KPATA" ||
str2 == "KPATK" ||
str2 == "KPAYC" ||
str2 == "KPEBE" ||
str2 == "KPEBH" ||
str2 == "KPEKT" ||
str2 == "KPEMA" ||
str2 == "KPEMC" ||
str2 == "KPECC" ||
str2 == "KPECT" ||
str2 == "KPEYC" ||
str2 == "KPY3A" ||
str2 == "KPY3E" ||
str2 == "KPYKC" ||
str2 == "KPYMM" ||
str2 == "KPYTA" ||
str2 == "KPYYC" ||
str2 == "KPY4A" ||
str2 == "KPY4E" ||
str2 == "KPY4Y" ||
str2 == "KYBAC" ||
str2 == "KYBBA" ||
str2 == "KYBE3" ||
str2 == "KY3EH" ||
str2 == "KYKAL" ||
str2 == "KYKAH" ||
str2 == "KYKLA" ||
str2 == "KYKCA" ||
str2 == "KYKTA" ||
str2 == "KYKYP" ||
str2 == "KYLAK" ||
str2 == "KYLAH" ||
str2 == "KYLAC" ||
str2 == "KYLAY" ||
str2 == "KYLA4" ||
str2 == "KYLEK" ||
str2 == "KYLEP" ||
str2 == "KYMAK" ||
str2 == "KYMA4" ||
str2 == "KYMKA" ||
str2 == "KYMCA" ||
str2 == "KYMYX" ||
str2 == "KYHAK" ||
str2 == "KYHAY" ||
str2 == "KYHEB" ||
str2 == "KYHKA" ||
str2 == "KYHCA" ||
str2 == "KYHTA" ||
str2 == "KYPAK" ||
str2 == "KYPAM" ||
str2 == "KYPAH" ||
str2 == "KYPAT" ||
str2 == "KYPBA" ||
str2 == "KYPEC" ||
str2 == "KYPKA" ||
str2 == "KYPMA" ||
str2 == "KYPHA" ||
str2 == "KYPCA" ||
str2 == "KYPCK" ||
str2 == "KYPTA" ||
str2 == "KYPYK" ||
str2 == "KYPYL" ||
str2 == "KYPYH" ||
str2 == "KYCTY" ||
str2 == "KYTAK" ||
str2 == "KYTAP" ||
str2 == "KYTAC" ||
str2 == "KYTBA" ||
str2 == "KYTEM" ||
str2 == "KYTEH" ||
str2 == "KYTEP" ||
str2 == "KYTPA" ||
str2 == "KYTY3" ||
str2 == "KYTYM" ||
str2 == "KYXBA" ||
str2 == "KYXTA" ||
str2 == "KY4AK" ||
str2 == "KY4AH" ||
str2 == "KY4AC" ||
str2 == "KY4EP" ||
str2 == "KY4KA" ||
str2 == "KY4MA" ||
str2 == "KY4YK" ||
str2 == "KY4YM" ||
str2 == "KXAMA" ||
str2 == "KXACA" ||
str2 == "KXATA" ||
str2 == "KXMEP" ||
str2 == "LAALE" ||
str2 == "LABAH" ||
str2 == "LABAC" ||
str2 == "LABKA" ||
str2 == "LABPA" ||
str2 == "LABYA" ||
str2 == "LABYP" ||
str2 == "LA3AP" ||
str2 == "LA3EP" ||
str2 == "LA3KA" ||
str2 == "LA3YH" ||
str2 == "LAKAH" ||
str2 == "LAKEH" ||
str2 == "LAMAX" ||
str2 == "LAMEP" ||
str2 == "LAMET" ||
str2 == "LAMEX" ||
str2 == "LAMYT" ||
str2 == "LAHKA" ||
str2 == "LAPEK" ||
str2 == "LAPPA" ||
str2 == "LAPCA" ||
str2 == "LACKA" ||
str2 == "LACTA" ||
str2 == "LACYH" ||
str2 == "LATAM" ||
str2 == "LATKA" ||
str2 == "LATTA" ||
str2 == "LATYK" ||
str2 == "LATYP" ||
str2 == "LAYMA" ||
str2 == "LAYME" ||
str2 == "LAYPA" ||
str2 == "LAXEC" ||
str2 == "LAXTA" ||
str2 == "LA4KA" ||
str2 == "LEBAK" ||
str2 == "LEBEE" ||
str2 == "LEBEK" ||
str2 == "LEBEH" ||
str2 == "LE3EH" ||
str2 == "LE3EP" ||
str2 == "LEKAH" ||
str2 == "LEKEH" ||
str2 == "LEKMA" ||
str2 == "LEKCA" ||
str2 == "LELEB" ||
str2 == "LELEK" ||
str2 == "LEMAH" ||
str2 == "LEMEP" ||
str2 == "LEMET" ||
str2 == "LEMEX" ||
str2 == "LEMKE" ||
str2 == "LEMMA" ||
str2 == "LEMME" ||
str2 == "LEMYP" ||
str2 == "LEHAY" ||
str2 == "LEHBA" ||
str2 == "LEHEH" ||
str2 == "LEHKA" ||
str2 == "LEHCK" ||
str2 == "LEHTA" ||
str2 == "LEPKA" ||
str2 == "LEPMA" ||
str2 == "LEPME" ||
str2 == "LEPHA" ||
str2 == "LEPYA" ||
str2 == "LEPXE" ||
str2 == "LECKA" ||
str2 == "LECHA" ||
str2 == "LECXA" ||
str2 == "LETAC" ||
str2 == "LETKA" ||
str2 == "LETTE" ||
str2 == "LETYH" ||
str2 == "LETYC" ||
str2 == "LEYKA" ||
str2 == "LEXEP" ||
str2 == "LE44E" ||
str2 == "LYAPA" ||
str2 == "LYBEH" ||
str2 == "LYBYA" ||
str2 == "LY3AH" ||
str2 == "LY3EH" ||
str2 == "LYKAH" ||
str2 == "LYKAC" ||
str2 == "LYKAY" ||
str2 == "LYKA4" ||
str2 == "LYKKA" ||
str2 == "LYLEA" ||
str2 == "LYLYA" ||
str2 == "LYMKA" ||
str2 == "LYHAP" ||
str2 == "LYHEK" ||
str2 == "LYHKA" ||
str2 == "LYHXA" ||
str2 == "LYCKA" ||
str2 == "LYCTA" ||
str2 == "LYXTA" ||
str2 == "LY4KA" ||
str2 == "LXACA" ||
str2 == "MAACA" ||
str2 == "MAAXA" ||
str2 == "MABKA" ||
str2 == "MAETA" ||
str2 == "MA3AP" ||
str2 == "MA3AC" ||
str2 == "MA3EP" ||
str2 == "MA3KA" ||
str2 == "MA3YP" ||
str2 == "MA3YT" ||
str2 == "MAKAM" ||
str2 == "MAKAP" ||
str2 == "MAKET" ||
str2 == "MAKKA" ||
str2 == "MAKKE" ||
str2 == "MAKCA" ||
str2 == "MAKTA" ||
str2 == "MAKYA" ||
str2 == "MALBA" ||
str2 == "MALEA" ||
str2 == "MALEK" ||
str2 == "MALEP" ||
str2 == "MALKA" ||
str2 == "MAMAH" ||
str2 == "MAMKA" ||
str2 == "MAMPE" ||
str2 == "MAHAC" ||
str2 == "MAHAT" ||
str2 == "MAHEP" ||
str2 == "MAHEC" ||
str2 == "MAHKA" ||
str2 == "MAHKE" ||
str2 == "MAHHA" ||
str2 == "MAHTA" ||
str2 == "MAHTY" ||
str2 == "MAHYL" ||
str2 == "MAHYC" ||
str2 == "MAH4A" ||
str2 == "MAPAL" ||
str2 == "MAPAH" ||
str2 == "MAPAT" ||
str2 == "MAPAX" ||
str2 == "MAPEK" ||
str2 == "MAPKA" ||
str2 == "MAPKE" ||
str2 == "MAPKC" ||
str2 == "MAPKY" ||
str2 == "MAPHA" ||
str2 == "MAPCA" ||
str2 == "MAPTA" ||
str2 == "MAPTE" ||
str2 == "MAPYH" ||
str2 == "MAPYX" ||
str2 == "MAPXA" ||
str2 == "MACAK" ||
str2 == "MACAH" ||
str2 == "MACET" ||
str2 == "MACKA" ||
str2 == "MACLA" ||
str2 == "MACCA" ||
str2 == "MACCE" ||
str2 == "MATAC" ||
str2 == "MATEB" ||
str2 == "MATKA" ||
str2 == "MATLE" ||
str2 == "MATPA" ||
str2 == "MATYA" ||
str2 == "MATXA" ||
str2 == "MAT4A" ||
str2 == "MAYKA" ||
str2 == "MAYLE" ||
str2 == "MAXAL" ||
str2 == "MAXAH" ||
str2 == "MAXAP" ||
str2 == "MAXA4" ||
str2 == "MAXMA" ||
str2 == "MAXPA" ||
str2 == "MAXTA" ||
str2 == "MA4BA" ||
str2 == "MA4TA" ||
str2 == "MBEPY" ||
str2 == "MEBAH" ||
str2 == "MEBAT" ||
str2 == "MEBEC" ||
str2 == "MEEBA" ||
str2 == "ME3AP" ||
str2 == "ME3EP" ||
str2 == "MEKAH" ||
str2 == "MEKKA" ||
str2 == "MELAK" ||
str2 == "MELAM" ||
str2 == "MELAH" ||
str2 == "MELAP" ||
str2 == "MELAC" ||
str2 == "MELEM" ||
str2 == "MELEH" ||
str2 == "MELEP" ||
str2 == "MELKA" ||
str2 == "MELY3" ||
str2 == "MEHAM" ||
str2 == "MEHAH" ||
str2 == "MEHEE" ||
str2 == "MEHEK" ||
str2 == "MEHEM" ||
str2 == "MEHEH" ||
str2 == "MEHKA" ||
str2 == "MEHKE" ||
str2 == "MEHCY" ||
str2 == "MEHTA" ||
str2 == "MEHYA" ||
str2 == "MEH4Y" ||
str2 == "MEPAH" ||
str2 == "MEPEK" ||
str2 == "MEPEH" ||
str2 == "MEPET" ||
str2 == "MEPKA" ||
str2 == "MEPKE" ||
str2 == "MEPKC" ||
str2 == "MEPLE" ||
str2 == "MEPHC" ||
str2 == "MEPPA" ||
str2 == "MEPTA" ||
str2 == "MECEH" ||
str2 == "MECCA" ||
str2 == "MECTA" ||
str2 == "MECTP" ||
str2 == "METAH" ||
str2 == "METKA" ||
str2 == "METLA" ||
str2 == "METPA" ||
str2 == "METYL" ||
str2 == "MET4E" ||
str2 == "MEXPA" ||
str2 == "ME4KA" ||
str2 == "ME4TA" ||
str2 == "MLABA" ||
str2 == "MLAKA" ||
str2 == "MLEEB" ||
str2 == "MHEMA" ||
str2 == "MPACA" ||
str2 == "MTAXA" ||
str2 == "MTE3A" ||
str2 == "MYABP" ||
str2 == "MYAPE" ||
str2 == "MY3EP" ||
str2 == "MYKAP" ||
str2 == "MYLAT" ||
str2 == "MYLEH" ||
str2 == "MYLLA" ||
str2 == "MYPAT" ||
str2 == "MYP3A" ||
str2 == "MYPKA" ||
str2 == "MYPCY" ||
str2 == "MYPYH" ||
str2 == "MYCAH" ||
str2 == "MYCAT" ||
str2 == "MYCCA" ||
str2 == "MYTEB" ||
str2 == "MYTEP" ||
str2 == "MYTYL" ||
str2 == "MYXA4" ||
str2 == "MY4KA" ||
str2 == "HABAL" ||
str2 == "HABAP" ||
str2 == "HABAT" ||
str2 == "HABEK" ||
str2 == "HABEL" ||
str2 == "HABEC" ||
str2 == "HABET" ||
str2 == "HABCE" ||
str2 == "HAECT" ||
str2 == "HA3EM" ||
str2 == "HAKA3" ||
str2 == "HAKAL" ||
str2 == "HAKAT" ||
str2 == "HAKPA" ||
str2 == "HALEH" ||
str2 == "HALET" ||
str2 == "HAMA3" ||
str2 == "HAMEK" ||
str2 == "HAMEL" ||
str2 == "HAMET" ||
str2 == "HAMHY" ||
str2 == "HAHAK" ||
str2 == "HAHKA" ||
str2 == "HAPBA" ||
str2 == "HAPEB" ||
str2 == "HAPE3" ||
str2 == "HAPEK" ||
str2 == "HAPCA" ||
str2 == "HAPTA" ||
str2 == "HAP4A" ||
str2 == "HACEK" ||
str2 == "HACEL" ||
str2 == "HACEP" ||
str2 == "HACCA" ||
str2 == "HACCE" ||
str2 == "HATAL" ||
str2 == "HATAH" ||
str2 == "HATEK" ||
str2 == "HATEC" ||
str2 == "HATPY" ||
str2 == "HATTA" ||
str2 == "HAYKA" ||
str2 == "HAYPA" ||
str2 == "HAYPY" ||
str2 == "HAXAL" ||
str2 == "HAXAP" ||
str2 == "HAXYP" ||
str2 == "HA4EB" ||
str2 == "HA4EL" ||
str2 == "HA4EC" ||
str2 == "HA4ET" ||
str2 == "HA4LA" ||
str2 == "HA4HY" ||
str2 == "HA4TY" ||
str2 == "HEAPX" ||
str2 == "HEBEP" ||
str2 == "HEBMA" ||
str2 == "HEBYC" ||
str2 == "HEE3E" ||
str2 == "HEELA" ||
str2 == "HEKA3" ||
str2 == "HEKAP" ||
str2 == "HEKEM" ||
str2 == "HEKMA" ||
str2 == "HEKCE" ||
str2 == "HELET" ||
str2 == "HEMAH" ||
str2 == "HEMET" ||
str2 == "HEMKA" ||
str2 == "HEHKA" ||
str2 == "HEPA3" ||
str2 == "HEPBA" ||
str2 == "HEPET" ||
str2 == "HEPKA" ||
str2 == "HEP4A" ||
str2 == "HECYH" ||
str2 == "HETEP" ||
str2 == "HETEC" ||
str2 == "HEYKA" ||
str2 == "HEYMA" ||
str2 == "HE4AC" ||
str2 == "HE4EM" ||
str2 == "HE4ET" ||
str2 == "HYAPE" ||
str2 == "HYKEP" ||
str2 == "HYKYC" ||
str2 == "HYMEA" ||
str2 == "HYMEH" ||
str2 == "HYMEP" ||
str2 == "HYPEH" ||
str2 == "HYPMA" ||
str2 == "HYTKA" ||
str2 == "HYTTA" ||
str2 == "PABEH" ||
str2 == "PAB3E" ||
str2 == "PAEKA" ||
str2 == "PA3BE" ||
str2 == "PA3EC" ||
str2 == "PA3YM" ||
str2 == "PAKAH" ||
str2 == "PAKKA" ||
str2 == "PAKCA" ||
str2 == "PAKYH" ||
str2 == "PAMAH" ||
str2 == "PAMEH" ||
str2 == "PAMKA" ||
str2 == "PAMLE" ||
str2 == "PAMYC" ||
str2 == "PAHAH" ||
str2 == "PAHEE" ||
str2 == "PAHEP" ||
str2 == "PAHET" ||
str2 == "PAHKA" ||
str2 == "PAHKE" ||
str2 == "PAHHE" ||
str2 == "PAHCE" ||
str2 == "PAHTE" ||
str2 == "PACCA" ||
str2 == "PACTP" ||
str2 == "PACTY" ||
str2 == "PACYL" ||
str2 == "PACXH" ||
str2 == "PATAH" ||
str2 == "PATAX" ||
str2 == "PATKA" ||
str2 == "PATKE" ||
str2 == "PAYCA" ||
str2 == "PAXAM" ||
str2 == "PAXEL" ||
str2 == "PA4EK" ||
str2 == "PEATC" ||
str2 == "PEBKA" ||
str2 == "PEBMA" ||
str2 == "PEBYH" ||
str2 == "PE3AK" ||
str2 == "PE3AH" ||
str2 == "PE3EK" ||
str2 == "PE3EH" ||
str2 == "PE3KA" ||
str2 == "PE3YC" ||
str2 == "PE34E" ||
str2 == "PEKKE" ||
str2 == "PEKTA" ||
str2 == "PEMAK" ||
str2 == "PEME3" ||
str2 == "PEMEP" ||
str2 == "PEMKE" ||
str2 == "PEHAH" ||
str2 == "PEHAP" ||
str2 == "PEHET" ||
str2 == "PEH3E" ||
str2 == "PEHHE" ||
str2 == "PEHTA" ||
str2 == "PEPYM" ||
str2 == "PECKA" ||
str2 == "PECCA" ||
str2 == "PECTY" ||
str2 == "PETEH" ||
str2 == "PETPA" ||
str2 == "PE4AE" ||
str2 == "PE4KA" ||
str2 == "PTY4Y" ||
str2 == "PY3MA" ||
str2 == "PYKAB" ||
str2 == "PYKAM" ||
str2 == "PYKET" ||
str2 == "PYLAH" ||
str2 == "PYLE3" ||
str2 == "PYLEK" ||
str2 == "PYLET" ||
str2 == "PYMEP" ||
str2 == "PYHET" ||
str2 == "PYHKA" ||
str2 == "PYCAK" ||
str2 == "PYCKA" ||
str2 == "PYCCE" ||
str2 == "PYTKA" ||
str2 == "PY4KA" ||
str2 == "CABAH" ||
str2 == "CABAP" ||
str2 == "CABAX" ||
str2 == "CABBA" ||
str2 == "CABEL" ||
str2 == "CABKA" ||
str2 == "CA3AH" ||
str2 == "CAKBA" ||
str2 == "CAKEH" ||
str2 == "CAKKC" ||
str2 == "CAKMA" ||
str2 == "CAKCE" ||
str2 == "CALAM" ||
str2 == "CALAC" ||
str2 == "CALAT" ||
str2 == "CALEM" ||
str2 == "CALEH" ||
str2 == "CALEX" ||
str2 == "CALKA" ||
str2 == "CALLE" ||
str2 == "CALMA" ||
str2 == "CALYH" ||
str2 == "CAMAH" ||
str2 == "CAMAP" ||
str2 == "CAMKA" ||
str2 == "CAMPA" ||
str2 == "CAMCA" ||
str2 == "CAMCE" ||
str2 == "CAMYA" ||
str2 == "CAMYM" ||
str2 == "CAMYP" ||
str2 == "CAMYC" ||
str2 == "CAHAA" ||
str2 == "CAHAT" ||
str2 == "CAHTA" ||
str2 == "CAPAH" ||
str2 == "CAPAY" ||
str2 == "CAPKA" ||
str2 == "CAPMA" ||
str2 == "CAPHA" ||
str2 == "CAPPA" ||
str2 == "CAPTA" ||
str2 == "CAPTP" ||
str2 == "CAPYM" ||
str2 == "CACCA" ||
str2 == "CACYH" ||
str2 == "CATAH" ||
str2 == "CATKA" ||
str2 == "CATPA" ||
str2 == "CAYLE" ||
str2 == "CAYHA" ||
str2 == "CAYPA" ||
str2 == "CAXAK" ||
str2 == "CAXAP" ||
str2 == "CAXTA" ||
str2 == "CBAHA" ||
str2 == "CBAPA" ||
str2 == "CBAPT" ||
str2 == "CBAXA" ||
str2 == "CBEHE" ||
str2 == "CBEPT" ||
str2 == "CBEPX" ||
str2 == "CBE4A" ||
str2 == "CBE4Y" ||
str2 == "CEAHC" ||
str2 == "CEAPA" ||
str2 == "CEBAK" ||
str2 == "CEBAH" ||
str2 == "CEBAP" ||
str2 == "CEBEP" ||
str2 == "CEBCK" ||
str2 == "CE3AM" ||
str2 == "CE3EP" ||
str2 == "CEKAM" ||
str2 == "CEKAC" ||
str2 == "CEKA4" ||
str2 == "CEKCT" ||
str2 == "CEKTA" ||
str2 == "CELAM" ||
str2 == "CELEH" ||
str2 == "CELEX" ||
str2 == "CELLA" ||
str2 == "CELMM" ||
str2 == "CEMAK" ||
str2 == "CEMAP" ||
str2 == "CEMEH" ||
str2 == "CEMME" ||
str2 == "CEHAK" ||
str2 == "CEHAT" ||
str2 == "CEHKT" ||
str2 == "CEHHA" ||
str2 == "CEHTA" ||
str2 == "CEPAM" ||
str2 == "CEPAH" ||
str2 == "CEPBE" ||
str2 == "CEPEE" ||
str2 == "CEPEH" ||
str2 == "CEPEP" ||
str2 == "CEPEC" ||
str2 == "CEPET" ||
str2 == "CEPKA" ||
str2 == "CEPLC" ||
str2 == "CEPHA" ||
str2 == "CEPPA" ||
str2 == "CEPPE" ||
str2 == "CEPYM" ||
str2 == "CECCE" ||
str2 == "CETAP" ||
str2 == "CETEP" ||
str2 == "CETKA" ||
str2 == "CETME" ||
str2 == "CEYH4" ||
str2 == "CEYTA" ||
str2 == "CEXEP" ||
str2 == "CE4KA" ||
str2 == "CKABP" ||
str2 == "CKALA" ||
str2 == "CKALY" ||
str2 == "CKAPA" ||
str2 == "CKAPH" ||
str2 == "CKAYT" ||
str2 == "CKA4Y" ||
str2 == "CKBEP" ||
str2 == "CKEHA" ||
str2 == "CKET4" ||
str2 == "CKLAL" ||
str2 == "CKPAL" ||
str2 == "CKPAM" ||
str2 == "CKYKA" ||
str2 == "CKYLA" ||
str2 == "CKYHC" ||
str2 == "CKYPA" ||
str2 == "CLABA" ||
str2 == "CLAKA" ||
str2 == "CLAMA" ||
str2 == "CLAHA" ||
str2 == "CLEBA" ||
str2 == "CLE3A" ||
str2 == "CLE4Y" ||
str2 == "CLYKA" ||
str2 == "CMALY" ||
str2 == "CMAPT" ||
str2 == "CMAXY" ||
str2 == "CMELA" ||
str2 == "CMEHA" ||
str2 == "CMEPA" ||
str2 == "CMEP4" ||
str2 == "CMETA" ||
str2 == "CMETC" ||
str2 == "CMETY" ||
str2 == "CME4Y" ||
str2 == "CMYTA" ||
str2 == "CMYYL" ||
str2 == "CMYXA" ||
str2 == "CHAHA" ||
str2 == "CHEEK" ||
str2 == "CHEKA" ||
str2 == "CHELL" ||
str2 == "CPA3Y" ||
str2 == "CCEKY" ||
str2 == "CTABP" ||
str2 == "CTA3A" ||
str2 == "CTAHC" ||
str2 == "CTAHY" ||
str2 == "CTAPA" ||
str2 == "CTAPE" ||
str2 == "CTAPK" ||
str2 == "CTAPP" ||
str2 == "CTAPT" ||
str2 == "CTEEH" ||
str2 == "CTEKA" ||
str2 == "CTEKY" ||
str2 == "CTELA" ||
str2 == "CTEHA" ||
str2 == "CTEPK" ||
str2 == "CTEPH" ||
str2 == "CTEPT" ||
str2 == "CTEPX" ||
str2 == "CTEXA" ||
str2 == "CTPA3" ||
str2 == "CTPAM" ||
str2 == "CTPAC" ||
str2 == "CTPAX" ||
str2 == "CTPEK" ||
str2 == "CTPYK" ||
str2 == "CTPYC" ||
str2 == "CTYPE" ||
str2 == "CYAPE" ||
str2 == "CYBAP" ||
str2 == "CYETA" ||
str2 == "CY3AB" ||
str2 == "CY3EM" ||
str2 == "CYKHA" ||
str2 == "CYKPE" ||
str2 == "CYLAK" ||
str2 == "CYLLA" ||
str2 == "CYLYK" ||
str2 == "CYMAK" ||
str2 == "CYMAX" ||
str2 == "CYMET" ||
str2 == "CYMKA" ||
str2 == "CYMMA" ||
str2 == "CYMYH" ||
str2 == "CYMYP" ||
str2 == "CYHHA" ||
str2 == "CYPAM" ||
str2 == "CYPAH" ||
str2 == "CYPAT" ||
str2 == "CYPLA" ||
str2 == "CYPMA" ||
str2 == "CYPHA" ||
str2 == "CYPPA" ||
str2 == "CYPCK" ||
str2 == "CYPTA" ||
str2 == "CYPYC" ||
str2 == "CYCAK" ||
str2 == "CYCEK" ||
str2 == "CYCYK" ||
str2 == "CYTPA" ||
str2 == "CYXYM" ||
str2 == "CY4AH" ||
str2 == "CY4KA" ||
str2 == "CXBAT" ||
str2 == "CXEMA" ||
str2 == "C4ETA" ||
str2 == "TABYH" ||
str2 == "TA3AK" ||
str2 == "TA3AH" ||
str2 == "TAKEB" ||
str2 == "TAKET" ||
str2 == "TAKKA" ||
str2 == "TAKLA" ||
str2 == "TAKHA" ||
str2 == "TAKCA" ||
str2 == "TAKTA" ||
str2 == "TALAK" ||
str2 == "TALAH" ||
str2 == "TALAP" ||
str2 == "TALAC" ||
str2 == "TALEB" ||
str2 == "TALEP" ||
str2 == "TALEC" ||
str2 == "TALCA" ||
str2 == "TAMAK" ||
str2 == "TAMAP" ||
str2 == "TAMKA" ||
str2 == "TAMYT" ||
str2 == "TAHAK" ||
str2 == "TAHAH" ||
str2 == "TAHKA" ||
str2 == "TAHHA" ||
str2 == "TAHTA" ||
str2 == "TAHYC" ||
str2 == "TAHXE" ||
str2 == "TAPA3" ||
str2 == "TAPAK" ||
str2 == "TAPAH" ||
str2 == "TAPAP" ||
str2 == "TAPAC" ||
str2 == "TAPKA" ||
str2 == "TAPLE" ||
str2 == "TAPHA" ||
str2 == "TAPTY" ||
str2 == "TACKA" ||
str2 == "TACMA" ||
str2 == "TACTY" ||
str2 == "TATKA" ||
str2 == "TATPA" ||
str2 == "TAYKE" ||
str2 == "TAYHC" ||
str2 == "TAXAH" ||
str2 == "TAXAY" ||
str2 == "TAXTA" ||
str2 == "TA4KA" ||
str2 == "TEATE" ||
str2 == "TEATP" ||
str2 == "TEBKP" ||
str2 == "TE3EK" ||
str2 == "TE3KA" ||
str2 == "TEKEC" ||
str2 == "TEKKE" ||
str2 == "TEKCT" ||
str2 == "TEKY4" ||
str2 == "TELAB" ||
str2 == "TELEM" ||
str2 == "TELKA" ||
str2 == "TELLE" ||
str2 == "TELYM" ||
str2 == "TELYP" ||
str2 == "TEMEH" ||
str2 == "TEM3A" ||
str2 == "TEMKA" ||
str2 == "TEMME" ||
str2 == "TEMHA" ||
str2 == "TEMHE" ||
str2 == "TEHAP" ||
str2 == "TEHEK" ||
str2 == "TEHET" ||
str2 == "TEH3A" ||
str2 == "TEHHY" ||
str2 == "TEPAC" ||
str2 == "TEPEK" ||
str2 == "TEPEM" ||
str2 == "TEPEH" ||
str2 == "TEPEC" ||
str2 == "TEPKA" ||
str2 == "TEPMA" ||
str2 == "TEPPA" ||
str2 == "TECAK" ||
str2 == "TECKA" ||
str2 == "TECLA" ||
str2 == "TECCE" ||
str2 == "TECTA" ||
str2 == "TETKA" ||
str2 == "TEXAC" ||
str2 == "TE4KA" ||
str2 == "T3YPA" ||
str2 == "TPAAT" ||
str2 == "TPABA" ||
str2 == "TPABE" ||
str2 == "TPAKT" ||
str2 == "TPAMM" ||
str2 == "TPAHE" ||
str2 == "TPAHC" ||
str2 == "TPACC" ||
str2 == "TPACT" ||
str2 == "TPATA" ||
str2 == "TPATP" ||
str2 == "TPAYH" ||
str2 == "TPAYP" ||
str2 == "TPA4Y" ||
str2 == "TPELA" ||
str2 == "TPEMA" ||
str2 == "TPEHK" ||
str2 == "TPEHT" ||
str2 == "TPEH4" ||
str2 == "TPECK" ||
str2 == "TPECC" ||
str2 == "TPECT" ||
str2 == "TPEYX" ||
str2 == "TPHKA" ||
str2 == "TPTCY" ||
str2 == "TPYHA" ||
str2 == "TPYCK" ||
str2 == "TPYXA" ||
str2 == "TCYKA" ||
str2 == "TYAPA" ||
str2 == "TY3LA" ||
str2 == "TYKAH" ||
str2 == "TYKA4" ||
str2 == "TYKEP" ||
str2 == "TYLEC" ||
str2 == "TYLY3" ||
str2 == "TYLYK" ||
str2 == "TYLYM" ||
str2 == "TYLYH" ||
str2 == "TYL4A" ||
str2 == "TYMAK" ||
str2 == "TYMAH" ||
str2 == "TYMAP" ||
str2 == "TYMAC" ||
str2 == "TYMEH" ||
str2 == "TYHKA" ||
str2 == "TYPAH" ||
str2 == "TYPAX" ||
str2 == "TYPA4" ||
str2 == "TYPEK" ||
str2 == "TYPEH" ||
str2 == "TYPKA" ||
str2 == "TYPKY" ||
str2 == "TYPLA" ||
str2 == "TYPMA" ||
str2 == "TYPHE" ||
str2 == "TYPPA" ||
str2 == "TYPPE" ||
str2 == "TYPYH" ||
str2 == "TYP4A" ||
str2 == "TYCCA" ||
str2 == "TYTAK" ||
str2 == "TYTEH" ||
str2 == "TYTTA" ||
str2 == "TY4EK" ||
str2 == "TY4KA" ||
str2 == "TXAHA" ||
str2 == "TXAHY" ||
str2 == "TXYAH" ||
str2 == "YAHKA" ||
str2 == "YAPTE" ||
str2 == "YBLEK" ||
str2 == "SCAMS" ||
str2 == "YBPAP" ||
str2 == "Y3AHC" ||
str2 == "Y3BAP" ||
str2 == "Y3EPK" ||
str2 == "Y3YPA" ||
str2 == "YKALA" ||
str2 == "YKA4Y" ||
str2 == "YKEPT" ||
str2 == "YKPAL" ||
str2 == "YKPYT" ||
str2 == "YKPYX" ||
str2 == "YKCYC" ||
str2 == "YKTXA" ||
str2 == "YLALA" ||
str2 == "YLE4Y" ||
str2 == "YLLAP" ||
str2 == "YLLAC" ||
str2 == "YME4Y" ||
str2 == "YHEEB" ||
str2 == "YHE4A" ||
str2 == "YH3EH" ||
str2 == "YHCET" ||
str2 == "YHTEP" ||
str2 == "YPABA" ||
str2 == "YPA3A" ||
str2 == "YPAKA" ||
str2 == "YPAHA" ||
str2 == "YPACA" ||
str2 == "YPBAH" ||
str2 == "YPEKE" ||
str2 == "YPEMA" ||
str2 == "YPKA4" ||
str2 == "YPLAH" ||
str2 == "YPMAH" ||
str2 == "YPCYL" ||
str2 == "YCAMA" ||
str2 == "YCEKY" ||
str2 == "YCLAP" ||
str2 == "YCMAH" ||
str2 == "YCTAB" ||
str2 == "YTEKY" ||
str2 == "YTEHA" ||
str2 == "YTEXA" ||
str2 == "YTPAM" ||
str2 == "YYPAC" ||
str2 == "YXBAL" ||
str2 == "YXBAT" ||
str2 == "Y4ETA" ||
str2 == "Y4KYP" ||
str2 == "XAA3E" ||
str2 == "XABEA" ||
str2 == "XA3AH" ||
str2 == "XA3AP" ||
str2 == "XAKAH" ||
str2 == "XAKAC" ||
str2 == "XAKEP" ||
str2 == "XAKKA" ||
str2 == "XALAL" ||
str2 == "XALAT" ||
str2 == "XALBA" ||
str2 == "XALEB" ||
str2 == "XALE3" ||
str2 == "XALLE" ||
str2 == "XALXA" ||
str2 == "XAMAH" ||
str2 == "XAMAC" ||
str2 == "XAM3A" ||
str2 == "XAMKA" ||
str2 == "XAMCA" ||
str2 == "XAMCE" ||
str2 == "XAMYH" ||
str2 == "XAHKA" ||
str2 == "XAHHA" ||
str2 == "XAHTA" ||
str2 == "XAHYM" ||
str2 == "XAPAL" ||
str2 == "XAPAP" ||
str2 == "XAPBA" ||
str2 == "XAPEC" ||
str2 == "XAP3A" ||
str2 == "XAPKE" ||
str2 == "XAPMC" ||
str2 == "XAPYH" ||
str2 == "XACAH" ||
str2 == "XACCA" ||
str2 == "XACCE" ||
str2 == "XATKA" ||
str2 == "XATMA" ||
str2 == "XAY3E" ||
str2 == "XAYCA" ||
str2 == "XAXAM" ||
str2 == "XAXMA" ||
str2 == "XA4EK" ||
str2 == "XBALA" ||
str2 == "XBA4Y" ||
str2 == "XEKEP" ||
str2 == "XEHHA" ||
str2 == "XEPEM" ||
str2 == "XEPEC" ||
str2 == "XECAH" ||
str2 == "XECCE" ||
str2 == "XEYPE" ||
str2 == "XLEMA" ||
str2 == "XLEHA" ||
str2 == "XLECT" ||
str2 == "XLYCA" ||
str2 == "XLYCT" ||
str2 == "XMAPA" ||
str2 == "XMYPA" ||
str2 == "XPECT" ||
str2 == "XPYCT" ||
str2 == "XYAHA" ||
str2 == "XYKAP" ||
str2 == "XYMCA" ||
str2 == "XYHHY" ||
str2 == "XYHTA" ||
str2 == "XYPAL" ||
str2 == "XYPMA" ||
str2 == "XYPTA" ||
str2 == "XYPYL" ||
str2 == "XYPXA" ||
str2 == "XYCTA" ||
str2 == "XYTAH" ||
str2 == "XYTPA" ||
str2 == "4ABEC" ||
str2 == "4AKAH" ||
str2 == "4AKBA" ||
str2 == "4AKPA" ||
str2 == "4ALKA" ||
str2 == "4ALLA" ||
str2 == "4ALMA" ||
str2 == "4AMEK" ||
str2 == "4AMPA" ||
str2 == "4AMYP" ||
str2 == "4AHAK" ||
str2 == "4AHAX" ||
str2 == "4AHKA" ||
str2 == "4AH4A" ||
str2 == "4APAX" ||
str2 == "4APEK" ||
str2 == "4APKA" ||
str2 == "4APL3" ||
str2 == "4ACAP" ||
str2 == "4ATAL" ||
str2 == "4ATAM" ||
str2 == "4ATEM" ||
str2 == "4AYLA" ||
str2 == "4AYPA" ||
str2 == "4AXAH" ||
str2 == "4EKAL" ||
str2 == "4EKAH" ||
str2 == "4EKTA" ||
str2 == "4ELAK" ||
str2 == "4ELEK" ||
str2 == "4ELEH" ||
str2 == "4ELKA" ||
str2 == "4ELMA" ||
str2 == "4ELHA" ||
str2 == "4EME3" ||
str2 == "4EMEP" ||
str2 == "4EMPA" ||
str2 == "4EMYP" ||
str2 == "4EHLA" ||
str2 == "4EHCY" ||
str2 == "4EPBA" ||
str2 == "4EPE3" ||
str2 == "4EPEK" ||
str2 == "4EPEH" ||
str2 == "4EPEC" ||
str2 == "4EPET" ||
str2 == "4EPEX" ||
str2 == "4EPTA" ||
str2 == "4EPTE" ||
str2 == "4EPTY" ||
str2 == "4EP4Y" ||
str2 == "4ECKA" ||
str2 == "4ECME" ||
str2 == "4ETKA" ||
str2 == "4ET4E" ||
str2 == "4E4AK" ||
str2 == "4E4EH" ||
str2 == "4E4ET" ||
str2 == "4LEHA" ||
str2 == "4MYKA" ||
str2 == "4YBAK" ||
str2 == "4YBAL" ||
str2 == "4YKAP" ||
str2 == "4YK4A" ||
str2 == "4YLAH" ||
str2 == "4YMAK" ||
str2 == "4YMAH" ||
str2 == "4YMKA" ||
str2 == "4YPAK" ||
str2 == "4YPAC" ||
str2 == "4YPEK" ||
str2 == "4YPKA" ||
str2 == "4YT4E" ||
str2 == "4YXHA" ||
str2 == "4Y4XE" ||
str2 == "KRSHN"
) {
	Krasisvo(passwordz, prizm_password, RSaddressPrizm);
}
if(str3 == "AAAA" || 
str3 == "BRA4" || 
str3 == "SCAM" || 
str3 == "HALF" || 
str3 == "ALYX" || 
str3 == "HALK" || 
str3 == "MYXA" || 
str3 == "CBET" || 
str3 == "A3AH" ||
str3 == "AKBA" ||
str3 == "ALAH" ||
str3 == "ALLA" ||
str3 == "AMAP" ||
str3 == "AHHA" ||
str3 == "AHYC" ||
str3 == "APEC" ||
str3 == "APKA" ||
str3 == "BAKA" ||
str3 == "BAHK" ||
str3 == "BEEP" ||
str3 == "BEHA" ||
str3 == "BEPA" ||
str3 == "BEPX" ||
str3 == "BKYC" ||
str3 == "BHYK" ||
str3 == "BPYH" ||
str3 == "BCEM" ||
str3 == "BCEX" ||
str3 == "3BYK" ||
str3 == "3EBA" ||
str3 == "3EBC" ||
str3 == "3HAK" ||
str3 == "KBAC" ||
str3 == "KEKC" ||
str3 == "KPEM" ||
str3 == "KPYT" ||
str3 == "KYPC" ||
str3 == "LAMA" ||
str3 == "LEHA" ||
str3 == "LEXA" ||
str3 == "LYKA" ||
str3 == "LYHA" ||
str3 == "MAKC" ||
str3 == "MAMA" ||
str3 == "MAPK" ||
str3 == "MAPC" ||
str3 == "MAPT" ||
str3 == "MATA" ||
str3 == "MEMC" ||
str3 == "MPAK" ||
str3 == "MYKA" ||
str3 == "MYXA" ||
str3 == "PAMA" ||
str3 == "PEKA" ||
str3 == "PEKC" ||
str3 == "PYKA" ||
str3 == "PYHA" ||
str3 == "PYCA" ||
str3 == "PYCE" ||
str3 == "PYCC" ||
str3 == "CBET" ||
str3 == "CEBA" ||
str3 == "CEKC" ||
str3 == "CEMA" ||
str3 == "CEPA" ||
str3 == "CCCP" ||
str3 == "CTAC" ||
str3 == "CTEH" ||
str3 == "CYKA" ||
str3 == "TAHK" ||
str3 == "TATY" ||
str3 == "TEMA" ||
str3 == "TPAX" ||
str3 == "TPEK" ||
str3 == "TPYC" ||
str3 == "TYPE" ||
str3 == "YKA3" ||
str3 == "YLET" ||
str3 == "YHAC" ||
str3 == "YPAL" ||
str3 == "YPAH" ||
str3 == "YTKA" ||
str3 == "XALK" ||
str3 == "XAXA" ||
str3 == "XLAM" ||
str3 == "XPAM" ||
str3 == "XPEH" ||
str3 == "HTML" ||
str3 == "XYAH" ||
str3 == "4EPT" ||
str3 == "4LEN" ||
str3 == "LUKA" || 
str3 == "MATA" || 
str3 == "MAMA" || 
str3 == "BEPX" || 
str3 == "BBBB" || 
str3 == "CCCC" || 
str3 == "DDDD" || 
str3 == "EEEE" || 
str3 == "FFFF" || 
str3 == "GGGG" || 
str3 == "HHHH" || 
str3 == "MMMM" || 
str3 == "NNNN" || 
str3 == "LLLL" || 
str3 == "FFFF" || 
str3 == "GGGG" || 
str3 == "HHHH" || 
str3 == "JJJJ" || 
str3 == "KKKK" || 
str3 == "LLLL" || 
str3 == "MMMM" || 
str3 == "NNNN" || 
str3 == "PPPP" || 
str3 == "QQQQ" || 
str3 == "RRRR" || 
str3 == "SSSS" || 
str3 == "TTTT" || 
str3 == "UUUU" || 
str3 == "VVVV" || 
str3 == "WWWW" || 
str3 == "XXXX" || 
str3 == "YYYY" || 
str3 == "ZZZZ" || 
str3 == "2222" || 
str3 == "3333" || 
str3 == "4444" || 
str3 == "5555" || 
str3 == "6666" || 
str3 == "7777" || 
str3 == "8888" || 
str3 == "9999") {
	Krasisvo(passwordz, prizm_password, RSaddressPrizm);
}

if( 
str2 == "SLAVA" ||
str2 == "PA3UM" ||
str2 == "BE4ER" ||
str2 == "STEAM" ||
str2 == "BLAGA" ||
str2 == "YDA4A" ||
str2 == "FAKES" ||
str2 == "FEMES" ||
str2 == "EDGES" ||
str2 == "DREAM" ||
str2 == "DUREX" ||
str2 == "DEEPS" ||
str2 == "DEBUG" ||
str2 == "CLANS" ||
str2 == "CHEAT" ||
str2 == "ABAMA" ||
str2 == "PUSSY" ||
str2 == "P7R4Q" ||
str2 == "Z8RJ5") {
		Krasisvo(passwordz, prizm_password, RSaddressPrizm);
}

} }



function GenerateRandomAdress(n, passwordz) {
	var fs = require('fs');
	for (i = 0; i < n; i++) {
       a = getRandomInt(1, 1000)
       prizm_password = gen_password(a);
	   PublicKeyPrizm = getPublicKeyPrizm(passwordz + ":" + prizm_password)
	   AccountId = getAccountId(PublicKeyPrizm)
	   RSaddressPrizm = getRSaddressPrizm(AccountId)
	   str2 = RSaddressPrizm.slice(21)
	   str3 = RSaddressPrizm.slice(6, 10)
	   Krasisvo(passwordz, prizm_password, RSaddressPrizm);
	}
	main();
} 



function main() {
	var readline = require('readline-sync');
	console.log("Генератор кошельков криптовалюты PRIZM");
	console.log("Разработчик - @Viktor_SMI - Telegram");
	console.log(" ");
	readline.question("Enter");
	console.log(" ");
	console.log("1. Сгенерировать N кошельков.");
	console.log("2. Рандомный красивый кошелек.");
	console.log("3. Сгенерировать красивый кошелек.");
	var a = readline.question("1/2/3 - ");
	if (a) {
	if (a == 1) {
console.log("Введите пароль:");
passwords = readline.question("");
console.log("Введите количество кошельков:");
n = readline.question("");
if (n != 0) {
GenerateRandomAdress(n, passwords);
} else {console.log("0 кошельков уже сгенерировано.")}
		}
	if (a == 2) { 
	console.log("Введите пароль:");
    passwords = readline.question("");
	GenerateRanAdress(passwords);
	}
	if (a == 3) {  
	getName();
	}
	} else {
        console.clear();
		console.log('Ошибка');
        console.log(' ');
		console.log(' ');
		console.log(' ');
		console.log(' ');
		console.log(' ');
		main();

} 
} 


main();