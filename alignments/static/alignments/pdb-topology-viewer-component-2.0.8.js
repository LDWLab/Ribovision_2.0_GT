function getCol(matrix, col){
    var column = [];
    for(var i=0; i<matrix.length; i++){
       column.push(matrix[i][col]);
    }
    return column;
 }


function hexToRgb(hex) {
    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
      r: parseInt(result[1], 16),
      g: parseInt(result[2], 16),
      b: parseInt(result[3], 16)
    } : null;
  }
  

  function hexToRgb1(hex) {
    var bigint = parseInt(hex, 16);
    var r = (bigint >> 16) & 255;
    var g = (bigint >> 8) & 255;
    var b = bigint & 255;

    return r + "," + g + "," + b;
}

function hexToRgb_GreenBlind(hex) {
    var bigint = parseInt(hex, 16);
    var r = (bigint >> 16) & 255;
    var g = (bigint >> 8) & 255;
    var b = bigint & 255;

    var r1 = Math.pow(r, 2.2);
    var g1 = Math.pow(g, 2.2);
    var b1 = Math.pow(b, 2.2);

    var R = Math.pow(0.02138 + 0.677 * g1 + 0.2802 * r1, 1 / 2.2);
    var B = Math.pow(0.02138 * (1 + g1 - r1) + 0.9572 * b1, 1 / 2.2);

    return R + "," + R + "," + B;
}

function hexToRgb_GreenBlind2(hex) {
    var bigint = parseInt(hex, 16);
    var r = (bigint >> 16) & 255;
    var g = (bigint >> 8) & 255;
    var b = bigint & 255;

    var r1 = Math.pow(r, 2.2);
    var g1 = Math.pow(g, 2.2);
    var b1 = Math.pow(b, 2.2);
    var R = Math.pow(0.003974 + 0.8806 * g1 + 0.1115 * r1, 1 / 2.2);
    var B = Math.pow(0.003974 * (1 - g1 + r1) + 0.9921 * b1, 1 / 2.2);


    return R + "," + R + "," + B;
}

function GreenBlind(r, g, b) {
    r = Math.pow(r, 2.2);
    g = Math.pow(g, 2.2);
    b = Math.pow(b, 2.2);
    var R = Math.pow(0.02138 + 0.677 * g + 0.2802 * r, 1 / 2.2);
    var B = Math.pow(0.02138 * (1 + g - r) + 0.9572 * b, 1 / 2.2);
    return [R, R, B];
  }

/**
 * pdb-topology-viewer
 * @version 2.0.0
 * @link https://github.com/PDBeurope/pdb-topology-viewer
 * @license Apache 2.0
 */
/**
 * Copyright 2020-2021 Mandar Deshpande <mandar@ebi.ac.uk>
 * European Bioinformatics Institute (EBI, http://www.ebi.ac.uk/)
 * European Molecular Biology Laboratory (EMBL, http://www.embl.de/)
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, 
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and 
 * limitations under the License.
 */


'use strict';

var palette = (function() {

  var proto = Array.prototype;
  var slice = function(arr, opt_begin, opt_end) {
    return proto.slice.apply(arr, proto.slice.call(arguments, 1));
  };

  var extend = function(arr, arr2) {
    return proto.push.apply(arr, arr2);
  };

  var function_type = typeof function() {};

  var INF = 1000000000;  // As far as we're concerned, that's infinity. ;)


  /**
   * Generate a colour palette from given scheme.
   *
   * If scheme argument is not a function it is passed to palettes.listSchemes
   * function (along with the number argument).  This may result in an array
   * of more than one available scheme.  If that is the case, scheme at
   * opt_index position is taken.
   *
   * This allows using different palettes for different data without having to
   * name the schemes specifically, for example:
   *
   *     palette_for_foo = palette('sequential', 10, 0);
   *     palette_for_bar = palette('sequential', 10, 1);
   *     palette_for_baz = palette('sequential', 10, 2);
   *
   * @param {!palette.SchemeType|string|palette.Palette} scheme Scheme to
   *     generate palette for.  Either a function constructed with
   *     palette.Scheme object, or anything that palette.listSchemes accepts
   *     as name argument.
   * @param {number} number Number of colours to return.  If negative, absolute
   *     value is taken and colours will be returned in reverse order.
   * @param {number=} opt_index If scheme is a name of a group or an array and
   *     results in more than one scheme, index of the scheme to use.  The
   *     index wraps around.
   * @param {...*} varargs Additional arguments to pass to palette or colour
   *     generator (if the chosen scheme uses those).
   * @return {Array<string>} Array of abs(number) 'RRGGBB' strings or null if
   *     no matching scheme was found.
   */
  var palette = function(scheme, number, opt_index, varargs) {
    number |= 0;
    if (number == 0) {
      return [];
    }

    if (typeof scheme !== function_type) {
      var arr = palette.listSchemes(
          /** @type {string|palette.Palette} */ (scheme), number);
      if (!arr.length) {
        return null;
      }
      scheme = arr[(opt_index || 0) % arr.length];
    }

    var args = slice(arguments, 2);
    args[0] = number;
    return scheme.apply(scheme, args);
  };


  /**
   * Returns a callable colour scheme object.
   *
   * Just after being created, the scheme has no colour palettes and no way of
   * generating any, thus generate method will return null.  To turn scheme
   * into a useful object, addPalette, addPalettes or setColorFunction methods
   * need to be used.
   *
   * To generate a colour palette with given number colours using function
   * returned by this method, just call it with desired number of colours.
   *
   * Since this function *returns* a callable object, it must *not* be used
   * with the new operator.
   *
   * @param {string} name Name of the scheme.
   * @param {string|!Array<string>=} opt_groups A group name or list of
   *     groups the scheme should be categorised under.  Three typical groups
   *     to use are 'qualitative', 'sequential' and 'diverging', but any
   *     groups may be created.
   * @return {!palette.SchemeType} A colour palette generator function, which
   *     in addition has methods and properties like a regular object.  Think
   *     of it as a callable object.
   */
  palette.Scheme = function(name, opt_groups) {
    /**
     * A map from a number to a colour palettes with given number of colours.
     * @type {!Object<number, palette.Palette>}
     */
    var palettes = {};

    /**
     * The biggest palette in palettes map.
     * @type {number}
     */
    var palettes_max = 0;

    /**
     * The smallest palette in palettes map.
     * @type {number}
     */
    var palettes_min = INF;

    var makeGenerator = function() {
      if (arguments.length <= 1) {
        return self.color_func.bind(self);
      } else {
        var args = slice(arguments);
        return function(x) {
          args[0] = x;
          return self.color_func.apply(self, args);
        };
      }
    };

    /**
     * Generate a colour palette from the scheme.
     *
     * If there was a palette added with addPalette (or addPalettes) with
     * enough colours, that palette will be used.  Otherwise, if colour
     * function has been set using setColorFunction method, that function will
     * be used to generate the palette.  Otherwise null is returned.
     *
     * @param {number} number Number of colours to return.  If negative,
     *     absolute value is taken and colours will be returned in reverse
     *     order.
     * @param {...*} varargs Additional arguments to pass to palette or colour
     *     generator (if the chosen scheme uses those).
     */
    var self = function(number, varargs) {
      number |= 0;
      if (!number) {
        return [];
      }

      var _number = number;
      number = Math.abs(number);

      if (number <= palettes_max) {
        for (var i = Math.max(number, palettes_min); !(i in palettes); ++i) {
          /* nop */
        }
        var colors = palettes[i];
        if (i > number) {
          var take_head =
              'shrinking_takes_head' in colors ?
              colors.shrinking_takes_head : self.shrinking_takes_head;
          if (take_head) {
            colors = colors.slice(0, number);
            i = number;
          } else {
            return palette.generate(
                function(x) { return colors[Math.round(x)]; },
                _number, 0, colors.length - 1);
          }
        }
        colors = colors.slice();
        if (_number < 0) {
          colors.reverse();
        }
        return colors;

      } else if (self.color_func) {
        return palette.generate(makeGenerator.apply(self, arguments),
                                _number, 0, 1, self.color_func_cyclic);

      } else {
        return null;
      }
    };

    /**
     * The name of the palette.
     * @type {string}
     */
    self.scheme_name = name;

    /**
     * A list of groups the palette belongs to.
     * @type {!Array<string>}
     */
    self.groups = opt_groups ?
      typeof opt_groups === 'string' ? [opt_groups] : opt_groups : [];

    /**
     * The biggest palette this scheme can generate.
     * @type {number}
     */
    self.max = 0;

    /**
     * The biggest palette this scheme can generate that is colour-blind
     * friendly.
     * @type {number}
     */
    self.cbf_max = INF;


    /**
     * Adds a colour palette to the colour scheme.
     *
     * @param {palette.Palette} palette An array of 'RRGGBB' strings
     *     representing the palette to add.
     * @param {boolean=} opt_is_cbf Whether the palette is colourblind friendly.
     */
    self.addPalette = function(palette, opt_is_cbf) {
      var len = palette.length;
      if (len) {
        palettes[len] = palette;
        palettes_min = Math.min(palettes_min, len);
        palettes_max = Math.max(palettes_max, len);
        self.max = Math.max(self.max, len);
        if (!opt_is_cbf && len != 1) {
          self.cbf_max = Math.min(self.cbf_max, len - 1);
        }
      }
    };

    /**
     * Adds number of colour palettes to the colour scheme.
     *
     * @param {palette.PalettesList} palettes A map or an array of colour
     *     palettes to add.  If map, i.e.  object, is used, properties should
     *     use integer property names.
     * @param {number=} opt_max Size of the biggest palette in palettes set.
     *     If not set, palettes must have a length property which will be used.
     * @param {number=} opt_cbf_max Size of the biggest palette which is still
     *     colourblind friendly.  1 by default.
     */
    self.addPalettes = function(palettes, opt_max, opt_cbf_max) {
      opt_max = opt_max || palettes.length;
      for (var i = 0; i < opt_max; ++i) {
        if (i in palettes) {
          self.addPalette(palettes[i], true);
        }
      }
      self.cbf_max = Math.min(self.cbf_max, opt_cbf_max || 1);
    };

    /**
     * Enable shrinking palettes taking head of the list of colours.
     *
     * When user requests n-colour palette but the smallest palette added with
     * addPalette (or addPalettes) is m-colour one (where n < m), n colours
     * across the palette will be returned.  For example:
     *     var ex = palette.Scheme('ex');
     *     ex.addPalette(['000000', 'bcbcbc', 'ffffff']);
     *     var pal = ex(2);
     *     // pal == ['000000', 'ffffff']
     *
     * This works for palettes where the distance between colours is
     * correlated to distance in the palette array, which is true in gradients
     * such as the one above.
     *
     * To turn this feature off shrinkByTakingHead can be set to true either
     * for all palettes in the scheme (if opt_idx is not given) or for palette
     * with given number of colours only.  In general, setting the option for
     * given palette overwrites whatever has been set for the scheme.  The
     * default, as described above, is false.
     *
     * Alternatively, the feature can be enabled by setting shrinking_takes_head
     * property for the palette Array or the scheme object.
     *
     * For example, all of the below give equivalent results:
     *     var pal = ['ff0000', '00ff00', '0000ff'];
     *
     *     var ex = palette.Scheme('ex');
     *     ex.addPalette(pal);               // ex(2) == ['ff0000', '0000ff']
     *     ex.shrinkByTakingHead(true);      // ex(2) == ['ff0000', '00ff00']
     *
     *     ex = palette.Scheme('ex');
     *     ex.addPalette(pal);               // ex(2) == ['ff0000', '0000ff']
     *     ex.shrinkByTakingHead(true, 3);   // ex(2) == ['ff0000', '00ff00']
     *
     *     ex = palette.Scheme('ex');
     *     ex.addPalette(pal);
     *     ex.addPalette(pal);               // ex(2) == ['ff0000', '0000ff']
     *     pal.shrinking_takes_head = true;  // ex(2) == ['ff0000', '00ff00']
     *
     * @param {boolean} enabled Whether to enable or disable the “shrinking
     *     takes head” feature.  It is disabled by default.
     * @param {number=} opt_idx If given, the “shrinking takes head” option
     *     for palette with given number of colours is set.  If such palette
     *     does not exist, nothing happens.
     */
    self.shrinkByTakingHead = function(enabled, opt_idx) {
      if (opt_idx !== void(0)) {
        if (opt_idx in palettes) {
          palettes[opt_idx].shrinking_takes_head = !!enabled;
        }
      } else {
        self.shrinking_takes_head = !!enabled;
      }
    };

    /**
     * Sets a colour generation function of the colour scheme.
     *
     * The function must accept a singe number argument whose value can be from
     * 0.0 to 1.0, and return a colour as an 'RRGGBB' string.  This function
     * will be used when generating palettes, i.e. if 11-colour palette is
     * requested, this function will be called with arguments 0.0, 0.1, …, 1.0.
     *
     * If the palette generated by the function is colourblind friendly,
     * opt_is_cbf should be set to true.
     *
     * In some cases, it is not desirable to reach 1.0 when generating
     * a palette.  This happens for hue-rainbows where the 0–1 range corresponds
     * to a 0°–360° range in hues, and since hue at 0° is the same as at 360°,
     * it's desired to stop short the end of the range when generating
     * a palette.  To accomplish this, opt_cyclic should be set to true.
     *
     * @param {palette.ColorFunction} func A colour generator function.
     * @param {boolean=} opt_is_cbf Whether palette generate with the function
     *     is colour-blind friendly.
     * @param {boolean=} opt_cyclic Whether colour at 0.0 is the same as the
     *     one at 1.0.
     */
    self.setColorFunction = function(func, opt_is_cbf, opt_cyclic) {
      self.color_func = func;
      self.color_func_cyclic = !!opt_cyclic;
      self.max = INF;
      if (!opt_is_cbf && self.cbf_max === INF) {
        self.cbf_max = 1;
      }
    };

    self.color = function(x, varargs) {
      if (self.color_func) {
        return self.color_func.apply(this, arguments);
      } else {
        return null;
      }
    };

    return self;
  };


  /**
   * Creates a new palette.Scheme and initialises it by calling addPalettes
   * method with the rest of the arguments.
   *
   * @param {string} name Name of the scheme.
   * @param {string|!Array<string>} groups A group name or list of group
   *     names the scheme belongs to.
   * @param {!Object<number, palette.Palette>|!Array<palette.Palette>}
   *     palettes A map or an array of colour palettes to add.  If map, i.e.
   *     object, is used, properties should use integer property names.
   * @param {number=} opt_max Size of the biggest palette in palettes set.
   *     If not set, palettes must have a length property which will be used.
   * @param {number=} opt_cbf_max Size of the biggest palette which is still
   *     colourblind friendly.  1 by default.
   * @return {!palette.SchemeType} A colour palette generator function, which
   *     in addition has methods and properties like a regular object.  Think
   *     of it as a callable object.
   */
  palette.Scheme.fromPalettes = function(name, groups,
                                         palettes, opt_max, opt_cbf_max) {
    var scheme = palette.Scheme(name, groups);
    scheme.addPalettes.apply(scheme, slice(arguments, 2));
    return scheme;
  };


  /**
   * Creates a new palette.Scheme and initialises it by calling
   * setColorFunction method with the rest of the arguments.
   *
   * @param {string} name Name of the scheme.
   * @param {string|!Array<string>} groups A group name or list of group
   *     names the scheme belongs to.
   * @param {palette.ColorFunction} func A colour generator function.
   * @param {boolean=} opt_is_cbf Whether palette generate with the function
   *     is colour-blind friendly.
   * @param {boolean=} opt_cyclic Whether colour at 0.0 is the same as the
   *     one at 1.0.
   * @return {!palette.SchemeType} A colour palette generator function, which
   *     in addition has methods and properties like a regular object.  Think
   *     of it as a callable object.
   */
  palette.Scheme.withColorFunction = function(name, groups,
                                              func, opt_is_cbf, opt_cyclic) {
    var scheme = palette.Scheme(name, groups);
    scheme.setColorFunction.apply(scheme, slice(arguments, 2));
    return scheme;
  };


  /**
   * A map of registered schemes.  Maps a scheme or group name to a list of
   * scheme objects.  Property name is either 'n-<name>' for single scheme
   * names or 'g-<name>' for scheme group names.
   *
   * @type {!Object<string, !Array<!Object>>}
   */
  var registered_schemes = {};


  /**
   * Registers a new colour scheme.
   *
   * @param {!palette.SchemeType} scheme The scheme to add.
   */
  palette.register = function(scheme) {
    registered_schemes['n-' + scheme.scheme_name] = [scheme];
    scheme.groups.forEach(function(g) {
      (registered_schemes['g-' + g] =
       registered_schemes['g-' + g] || []).push(scheme);
    });
    (registered_schemes['g-all'] =
       registered_schemes['g-all'] || []).push(scheme);
  };


  /**
   * List all schemes that match given name and number of colours.
   *
   * name argument can be either a string or an array of strings.  In the
   * former case, the function acts as if the argument was an array with name
   * as a single argument (i.e. “palette.listSchemes('foo')” is exactly the same
   * as “palette.listSchemes(['foo'])”).
   *
   * Each name can be either name of a palette (e.g. 'tol-sq' for Paul Tol's
   * sequential palette), or a name of a group (e.g. 'sequential' for all
   * sequential palettes).  Name can therefore map to a single scheme or
   * several schemes.
   *
   * Furthermore, name can be suffixed with '-cbf' to indicate that only
   * schemes that are colourblind friendly should be returned.  For example,
   * 'rainbow' returns a HSV rainbow scheme, but because it is not colourblind
   * friendly, 'rainbow-cbf' returns no schemes.
   *
   * Some schemes may produce colourblind friendly palettes for some number of
   * colours.  For example ColorBrewer's Dark2 scheme is colourblind friendly
   * if no more than 3 colours are generated.  If opt_number is not specified,
   * 'qualitative-cbf' will include 'cb-Dark2' but if opt_number is given as,
   * say, 5 it won't.
   *
   * Name can also be 'all' which will return all registered schemes.
   * Naturally, 'all-cbf' will return all colourblind friendly schemes.
   *
   * Schemes are added to the library using palette.register.  Schemes are
   * created using palette.Scheme function.  By default, the following schemes
   * are available:
   *
   *     Name            Description
   *     --------------  -----------------------------------------------------
   *     tol             Paul Tol's qualitative scheme, cbf, max 12 colours.
   *     tol-dv          Paul Tol's diverging scheme, cbf.
   *     tol-sq          Paul Tol's sequential scheme, cbf.
   *     tol-rainbow     Paul Tol's qualitative scheme, cbf.
   *
   *     rainbow         A rainbow palette.
   *
   *     cb-YlGn         ColorBrewer sequential schemes.
   *     cb-YlGnBu
   *     cb-GnBu
   *     cb-BuGn
   *     cb-PuBuGn
   *     cb-PuBu
   *     cb-BuPu
   *     cb-RdPu
   *     cb-PuRd
   *     cb-OrRd
   *     cb-YlOrRd
   *     cb-YlOrBr
   *     cb-Purples
   *     cb-Blues
   *     cb-Greens
   *     cb-Oranges
   *     cb-Reds
   *     cb-Greys
   *
   *     cb-PuOr         ColorBrewer diverging schemes.
   *     cb-BrBG
   *     cb-PRGn
   *     cb-PiYG
   *     cb-RdBu
   *     cb-RdGy
   *     cb-RdYlBu
   *     cb-Spectral
   *     cb-RdYlGn
   *
   *     cb-Accent       ColorBrewer qualitative schemes.
   *     cb-Dark2
   *     cb-Paired
   *     cb-Pastel1
   *     cb-Pastel2
   *     cb-Set1
   *     cb-Set2
   *     cb-Set3
   *
   *     sol-base        Solarized base colours.
   *     sol-accent      Solarized accent colours.
   *
   * The following groups are also available by default:
   *
   *     Name            Description
   *     --------------  -----------------------------------------------------
   *     all             All registered schemes.
   *     sequential      All sequential schemes.
   *     diverging       All diverging schemes.
   *     qualitative     All qualitative schemes.
   *     cb-sequential   All ColorBrewer sequential schemes.
   *     cb-diverging    All ColorBrewer diverging schemes.
   *     cb-qualitative  All ColorBrewer qualitative schemes.
   *
   * You can read more about Paul Tol's palettes at http://www.sron.nl/~pault/.
   * You can read more about ColorBrewer at http://colorbrewer2.org.
   *
   * @param {string|!Array<string>} name A name of a colour scheme, of
   *     a group of colour schemes, or an array of any of those.
   * @param {number=} opt_number When requesting only colourblind friendly
   *     schemes, number of colours the scheme must provide generating such
   *     that the palette is still colourblind friendly.  2 by default.
   * @return {!Array<!palette.SchemeType>} An array of colour scheme objects
   *     matching the criteria.  Sorted by scheme name.
   */
  palette.listSchemes = function(name, opt_number) {
    if (!opt_number) {
      opt_number = 2;
    } else if (opt_number < 0) {
      opt_number = -opt_number;
    }

    var ret = [];
    (typeof name === 'string' ? [name] : name).forEach(function(n) {
      var cbf = n.substring(n.length - 4) === '-cbf';
      if (cbf) {
        n = n.substring(0, n.length - 4);
      }
      var schemes =
          registered_schemes['g-' + n] ||
          registered_schemes['n-' + n] ||
          [];
      for (var i = 0, scheme; (scheme = schemes[i]); ++i) {
        if ((cbf ? scheme.cbf : scheme.max) >= opt_number) {
          ret.push(scheme);
        }
      }
    });

    ret.sort(function(a, b) {
      return a.scheme_name >= b.scheme_name ?
        a.scheme_name > b.scheme_name ? 1 : 0 : -1;
    });
    return ret;
  };


  /**
   * Generates a palette using given colour generating function.
   *
   * The color_func callback must accept a singe number argument whose value
   * can vary from 0.0 to 1.0 (or in general from opt_start to opt_end), and
   * return a colour as an 'RRGGBB' string.  This function will be used when
   * generating palettes, i.e. if 11-colour palette is requested, this
   * function will be called with arguments 0.0, 0.1, …, 1.0.
   *
   * In some cases, it is not desirable to reach 1.0 when generating
   * a palette.  This happens for hue-rainbows where the 0–1 range corresponds
   * to a 0°–360° range in hues, and since hue at 0° is the same as at 360°,
   * it's desired to stop short the end of the range when generating
   * a palette.  To accomplish this, opt_cyclic should be set to true.
   *
   * opt_start and opt_end may be used to change the range the colour
   * generation function is called with.  opt_end may be less than opt_start
   * which will case to traverse the range in reverse.  Another way to reverse
   * the palette is requesting negative number of colours.  The two methods do
   * not always lead to the same results (especially if opt_cyclic is set).
   *
   * @param {palette.ColorFunction} color_func A colours generating callback
   *     function.
   * @param {number} number Number of colours to generate in the palette.  If
   *     number is negative, colours in the palette will be reversed.  If only
   *     one colour is requested, colour at opt_start will be returned.
   * @param {number=} opt_start Optional starting point for the palette
   *     generation function.  Zero by default.
   * @param {number=} opt_end Optional ending point for the palette generation
   *     function.  One by default.
   * @param {boolean=} opt_cyclic If true, function will assume colour at
   *     point opt_start is the same as one at opt_end.
   * @return {palette.Palette} An array of 'RRGGBB' colours.
   */
  palette.generate = function(color_func, number, opt_start, opt_end,
                              opt_cyclic) {
    if (Math.abs(number) < 1) {
      return [];
    }

    opt_start = opt_start === void(0) ? 0 : opt_start;
    opt_end = opt_end === void(0) ? 1 : opt_end;

    if (Math.abs(number) < 2) {
      return [color_func(opt_start)];
    }

    var i = Math.abs(number);
    var x = opt_start;
    var ret = [];
    var step = (opt_end - opt_start) / (opt_cyclic ? i : (i - 1));

    for (; --i >= 0; x += step) {
      ret.push(color_func(x));
    }
    if (number < 0) {
      ret.reverse();
    }
    return ret;
  };


  /**
   * Clamps value to [0, 1] range.
   * @param {number} v Number to limit value of.
   * @return {number} If v is inside of [0, 1] range returns v, otherwise
   *     returns 0 or 1 depending which side of the range v is closer to.
   */
  var clamp = function(v) {
    return v > 0 ? (v < 1 ? v : 1) : 0;
  };

  /**
   * Converts r, g, b triple into RRGGBB hex representation.
   * @param {number} r Red value of the colour in the range [0, 1].
   * @param {number} g Green value of the colour in the range [0, 1].
   * @param {number} b Blue value of the colour in the range [0, 1].
   * @return {string} A lower-case RRGGBB representation of the colour.
   */
  palette.rgbColor = function(r, g, b) {
    return [r, g, b].map(function(v) {
      v = Number(Math.round(clamp(v) * 255)).toString(16);
      return v.length == 1 ? '0' + v : v;
    }).join('');
  };

  /**
   * Converts a linear r, g, b triple into RRGGBB hex representation.
   * @param {number} r Linear red value of the colour in the range [0, 1].
   * @param {number} g Linear green value of the colour in the range [0, 1].
   * @param {number} b Linear blue value of the colour in the range [0, 1].
   * @return {string} A lower-case RRGGBB representation of the colour.
   */
  palette.linearRgbColor = function(r, g, b) {
    // http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_RGB.html
    return [r, g, b].map(function(v) {
      v = clamp(v);
      if (v <= 0.0031308) {
        v = 12.92 * v;
      } else {
        v = 1.055 * Math.pow(v, 1 / 2.4) - 0.055;
      }
      v = Number(Math.round(v * 255)).toString(16);
      return v.length == 1 ? '0' + v : v;
    }).join('');
  };

  /**
   * Converts an HSV colours to RRGGBB hex representation.
   * @param {number} h Hue in the range [0, 1].
   * @param {number=} opt_s Saturation in the range [0, 1].  One by default.
   * @param {number=} opt_v Value in the range [0, 1].  One by default.
   * @return {string} An RRGGBB representation of the colour.
   */
  palette.hsvColor = function(h, opt_s, opt_v) {
    h *= 6;
    var s = opt_s === void(0) ? 1 : clamp(opt_s);
    var v = opt_v === void(0) ? 1 : clamp(opt_v);
    var x = v * (1 - s * Math.abs(h % 2 - 1));
    var m = v * (1 - s);
    switch (Math.floor(h) % 6) {
    case 0: return palette.rgbColor(v, x, m);
    case 1: return palette.rgbColor(x, v, m);
    case 2: return palette.rgbColor(m, v, x);
    case 3: return palette.rgbColor(m, x, v);
    case 4: return palette.rgbColor(x, m, v);
    default: return palette.rgbColor(v, m, x);
    }
  };

  palette.register(palette.Scheme.withColorFunction(
    'rainbow', 'qualitative', palette.hsvColor, false, true));

  return palette;
})();

/** @typedef {function(number): string} */
palette.ColorFunction;

/** @typedef {!Array<string>} */
palette.Palette;

/** @typedef {!Object<number, palette.Palette>|!Array<palette.Palette>} */
palette.PalettesList;

/**
 * @typedef {
 *   function(number, ...?): Array<string>|
 *   {
 *     scheme_name: string,
 *     groups: !Array<string>,
 *     max: number,
 *     cbf_max: number,
 *     addPalette: function(!Array<string>, boolean=),
 *     addPalettes: function(palette.PalettesList, number=, number=),
 *     shrinkByTakingHead: function(boolean, number=),
 *     setColorFunction: function(palette.ColorFunction, boolean=, boolean=),
 *     color: function(number, ...?): ?string}}
 */
palette.SchemeType;


/* mpn65 palette start here. ************************************************/

/* The ‘mpn65’ palette is designed for systems which show many graphs which
   don’t have custom colour palettes chosen by humans or if number of necessary
   colours isn’t know a priori. */

(function() {
  var scheme = palette.Scheme.fromPalettes('mpn65', 'qualitative', [[
    'ff0029', '377eb8', '66a61e', '984ea3', '00d2d5', 'ff7f00', 'af8d00',
    '7f80cd', 'b3e900', 'c42e60', 'a65628', 'f781bf', '8dd3c7', 'bebada',
    'fb8072', '80b1d3', 'fdb462', 'fccde5', 'bc80bd', 'ffed6f', 'c4eaff',
    'cf8c00', '1b9e77', 'd95f02', 'e7298a', 'e6ab02', 'a6761d', '0097ff',
    '00d067', '000000', '252525', '525252', '737373', '969696', 'bdbdbd',
    'f43600', '4ba93b', '5779bb', '927acc', '97ee3f', 'bf3947', '9f5b00',
    'f48758', '8caed6', 'f2b94f', 'eff26e', 'e43872', 'd9b100', '9d7a00',
    '698cff', 'd9d9d9', '00d27e', 'd06800', '009f82', 'c49200', 'cbe8ff',
    'fecddf', 'c27eb6', '8cd2ce', 'c4b8d9', 'f883b0', 'a49100', 'f48800',
    '27d0df', 'a04a9b'
  ]]);
  scheme.shrinkByTakingHead(true);
  palette.register(scheme);
})();

/* Paul Tol's schemes start here. *******************************************/
/* See http://www.sron.nl/~pault/ */

(function() {
  var rgb = palette.rgbColor;

  /**
   * Calculates value of a polynomial at given point.
   * For example, poly(x, 1, 2, 3) calculates value of “1 + 2*x + 3*X²”.
   * @param {number} x Value to calculate polynomial for.
   * @param {...number} varargs Coefficients of the polynomial specified in
   *     the order of rising powers of x including constant as the first
   *     variable argument.
   */
  var poly = function(x, varargs) {
    var i = arguments.length - 1, n = arguments[i];
    while (i > 1) {
      n = n * x + arguments[--i];
    }
    return n;
  };

  /**
   * Calculate approximate value of error function with maximum error of 0.0005.
   * See <https://en.wikipedia.org/wiki/Error_function>.
   * @param {number} x Argument of the error function.
   * @return {number} Value of error function for x.
   */
  var erf = function(x) {
    // https://en.wikipedia.org/wiki/Error_function#Approximation_with_elementary_functions
    // This produces a maximum error of 0.0005 which is more then we need.  In
    // the worst case, that error is multiplied by four and then divided by two
    // before being multiplied by 255, so in the end, the error is multiplied by
    // 510 which produces 0.255 which is less than a single colour step.
    var y = poly(Math.abs(x), 1, 0.278393, 0.230389, 0.000972, 0.078108);
    y *= y; // y^2
    y *= y; // y^4
    y = 1 - 1 / y;
    return x < 0 ? -y : y;
  };

  palette.register(palette.Scheme.fromPalettes('tol', 'qualitative', [
    ['4477aa'],
    ['4477aa', 'cc6677'],
    ['4477aa', 'ddcc77', 'cc6677'],
    ['4477aa', '117733', 'ddcc77', 'cc6677'],
    ['332288', '88ccee', '117733', 'ddcc77', 'cc6677'],
    ['332288', '88ccee', '117733', 'ddcc77', 'cc6677', 'aa4499'],
    ['332288', '88ccee', '44aa99', '117733', 'ddcc77', 'cc6677', 'aa4499'],
    ['332288', '88ccee', '44aa99', '117733', '999933', 'ddcc77', 'cc6677',
     'aa4499'],
    ['332288', '88ccee', '44aa99', '117733', '999933', 'ddcc77', 'cc6677',
     '882255', 'aa4499'],
    ['332288', '88ccee', '44aa99', '117733', '999933', 'ddcc77', '661100',
     'cc6677', '882255', 'aa4499'],
    ['332288', '6699cc', '88ccee', '44aa99', '117733', '999933', 'ddcc77',
     '661100', 'cc6677', '882255', 'aa4499'],
    ['332288', '6699cc', '88ccee', '44aa99', '117733', '999933', 'ddcc77',
     '661100', 'cc6677', 'aa4466', '882255', 'aa4499']
  ], 12, 12));

  /**
   * Calculates a colour along Paul Tol's sequential colours axis.
   * See <http://www.sron.nl/~pault/colourschemes.pdf> figure 7 and equation 1.
   * @param {number} x Position of the colour on the axis in the [0, 1] range.
   * @return {string} An RRGGBB representation of the colour.
   */
  palette.tolSequentialColor = function(x) {
    return rgb(1 - 0.392 * (1 + erf((x - 0.869) / 0.255)),
               1.021 - 0.456 * (1 + erf((x - 0.527) / 0.376)),
               1 - 0.493 * (1 + erf((x - 0.272) / 0.309)));
  };

  palette.register(palette.Scheme.withColorFunction(
    'tol-sq', 'sequential', palette.tolSequentialColor, true));

  /**
   * Calculates a colour along Paul Tol's diverging colours axis.
   * See <http://www.sron.nl/~pault/colourschemes.pdf> figure 8 and equation 2.
   * @param {number} x Position of the colour on the axis in the [0, 1] range.
   * @return {string} An RRGGBB representation of the colour.
   */
  palette.tolDivergingColor = function(x) {
    var g = poly(x, 0.572, 1.524, -1.811) / poly(x, 1, -0.291, 0.1574);
    return rgb(poly(x, 0.235, -2.13, 26.92, -65.5, 63.5, -22.36),
               g * g,
               1 / poly(x, 1.579, -4.03, 12.92, -31.4, 48.6, -23.36));
  };

  palette.register(palette.Scheme.withColorFunction(
    'tol-dv', 'diverging', palette.tolDivergingColor, true));

  /**
   * Calculates a colour along Paul Tol's rainbow colours axis.
   * See <http://www.sron.nl/~pault/colourschemes.pdf> figure 13 and equation 3.
   * @param {number} x Position of the colour on the axis in the [0, 1] range.
   * @return {string} An RRGGBB representation of the colour.
   */
  palette.tolRainbowColor = function(x) {
    return rgb(poly(x, 0.472, -0.567, 4.05) / poly(x, 1, 8.72, -19.17, 14.1),
               poly(x, 0.108932, -1.22635, 27.284, -98.577, 163.3, -131.395,
                    40.634),
               1 / poly(x, 1.97, 3.54, -68.5, 243, -297, 125));
  };

  palette.register(palette.Scheme.withColorFunction(
    'tol-rainbow', 'qualitative', palette.tolRainbowColor, true));
})();

(function() {
    /*
     * Those are not really designed to be used in graphs, but we're keeping
     * them here in case someone cares.
     */
    palette.register(palette.Scheme.fromPalettes('sol-base', 'sequential', [
      ['002b36', '073642', '586e75', '657b83', '839496', '93a1a1', 'eee8d5',
       'fdf6e3']
    ], 1, 8));
    palette.register(palette.Scheme.fromPalettes('sol-accent', 'qualitative', [
      ['b58900', 'cb4b16', 'dc322f', 'd33682', '6c71c4', '268bd2', '2aa198',
       '859900']
    ]));
  })();
  
  
  /* ColorBrewer colour schemes start here. ***********************************/
  /* See http://colorbrewer2.org/ */


(function() {
  var schemes = {
    YlGn: {
      type: 'sequential',
      cbf: 42,
      3: ['f7fcb9', 'addd8e', '31a354'],
      4: ['ffffcc', 'c2e699', '78c679', '238443'],
      5: ['ffffcc', 'c2e699', '78c679', '31a354', '006837'],
      6: ['ffffcc', 'd9f0a3', 'addd8e', '78c679', '31a354', '006837'],
      7: ['ffffcc', 'd9f0a3', 'addd8e', '78c679', '41ab5d', '238443',
          '005a32'],
      8: ['ffffe5', 'f7fcb9', 'd9f0a3', 'addd8e', '78c679', '41ab5d',
          '238443', '005a32'],
      9: ['ffffe5', 'f7fcb9', 'd9f0a3', 'addd8e', '78c679', '41ab5d',
          '238443', '006837', '004529']
    },
    YlGnBu: {
      type: 'sequential',
      cbf: 42,
      3: ['edf8b1', '7fcdbb', '2c7fb8'],
      4: ['ffffcc', 'a1dab4', '41b6c4', '225ea8'],
      5: ['ffffcc', 'a1dab4', '41b6c4', '2c7fb8', '253494'],
      6: ['ffffcc', 'c7e9b4', '7fcdbb', '41b6c4', '2c7fb8', '253494'],
      7: ['ffffcc', 'c7e9b4', '7fcdbb', '41b6c4', '1d91c0', '225ea8',
          '0c2c84'],
      8: ['ffffd9', 'edf8b1', 'c7e9b4', '7fcdbb', '41b6c4', '1d91c0',
          '225ea8', '0c2c84'],
      9: ['ffffd9', 'edf8b1', 'c7e9b4', '7fcdbb', '41b6c4', '1d91c0',
          '225ea8', '253494', '081d58']
    },
    GnBu: {
      type: 'sequential',
      cbf: 42,
      3: ['e0f3db', 'a8ddb5', '43a2ca'],
      4: ['f0f9e8', 'bae4bc', '7bccc4', '2b8cbe'],
      5: ['f0f9e8', 'bae4bc', '7bccc4', '43a2ca', '0868ac'],
      6: ['f0f9e8', 'ccebc5', 'a8ddb5', '7bccc4', '43a2ca', '0868ac'],
      7: ['f0f9e8', 'ccebc5', 'a8ddb5', '7bccc4', '4eb3d3', '2b8cbe',
          '08589e'],
      8: ['f7fcf0', 'e0f3db', 'ccebc5', 'a8ddb5', '7bccc4', '4eb3d3',
          '2b8cbe', '08589e'],
      9: ['f7fcf0', 'e0f3db', 'ccebc5', 'a8ddb5', '7bccc4', '4eb3d3',
          '2b8cbe', '0868ac', '084081']
    },
    BuGn: {
      type: 'sequential',
      cbf: 42,
      3: ['e5f5f9', '99d8c9', '2ca25f'],
      4: ['edf8fb', 'b2e2e2', '66c2a4', '238b45'],
      5: ['edf8fb', 'b2e2e2', '66c2a4', '2ca25f', '006d2c'],
      6: ['edf8fb', 'ccece6', '99d8c9', '66c2a4', '2ca25f', '006d2c'],
      7: ['edf8fb', 'ccece6', '99d8c9', '66c2a4', '41ae76', '238b45',
          '005824'],
      8: ['f7fcfd', 'e5f5f9', 'ccece6', '99d8c9', '66c2a4', '41ae76',
          '238b45', '005824'],
      9: ['f7fcfd', 'e5f5f9', 'ccece6', '99d8c9', '66c2a4', '41ae76',
          '238b45', '006d2c', '00441b']
    },
    PuBuGn: {
      type: 'sequential',
      cbf: 42,
      3: ['ece2f0', 'a6bddb', '1c9099'],
      4: ['f6eff7', 'bdc9e1', '67a9cf', '02818a'],
      5: ['f6eff7', 'bdc9e1', '67a9cf', '1c9099', '016c59'],
      6: ['f6eff7', 'd0d1e6', 'a6bddb', '67a9cf', '1c9099', '016c59'],
      7: ['f6eff7', 'd0d1e6', 'a6bddb', '67a9cf', '3690c0', '02818a',
          '016450'],
      8: ['fff7fb', 'ece2f0', 'd0d1e6', 'a6bddb', '67a9cf', '3690c0',
          '02818a', '016450'],
      9: ['fff7fb', 'ece2f0', 'd0d1e6', 'a6bddb', '67a9cf', '3690c0',
          '02818a', '016c59', '014636']
    },
    PuBu: {
      type: 'sequential',
      cbf: 42,
      3: ['ece7f2', 'a6bddb', '2b8cbe'],
      4: ['f1eef6', 'bdc9e1', '74a9cf', '0570b0'],
      5: ['f1eef6', 'bdc9e1', '74a9cf', '2b8cbe', '045a8d'],
      6: ['f1eef6', 'd0d1e6', 'a6bddb', '74a9cf', '2b8cbe', '045a8d'],
      7: ['f1eef6', 'd0d1e6', 'a6bddb', '74a9cf', '3690c0', '0570b0',
          '034e7b'],
      8: ['fff7fb', 'ece7f2', 'd0d1e6', 'a6bddb', '74a9cf', '3690c0',
          '0570b0', '034e7b'],
      9: ['fff7fb', 'ece7f2', 'd0d1e6', 'a6bddb', '74a9cf', '3690c0',
          '0570b0', '045a8d', '023858']
    },
    BuPu: {
      type: 'sequential',
      cbf: 42,
      3: ['e0ecf4', '9ebcda', '8856a7'],
      4: ['edf8fb', 'b3cde3', '8c96c6', '88419d'],
      5: ['edf8fb', 'b3cde3', '8c96c6', '8856a7', '810f7c'],
      6: ['edf8fb', 'bfd3e6', '9ebcda', '8c96c6', '8856a7', '810f7c'],
      7: ['edf8fb', 'bfd3e6', '9ebcda', '8c96c6', '8c6bb1', '88419d',
          '6e016b'],
      8: ['f7fcfd', 'e0ecf4', 'bfd3e6', '9ebcda', '8c96c6', '8c6bb1',
          '88419d', '6e016b'],
      9: ['f7fcfd', 'e0ecf4', 'bfd3e6', '9ebcda', '8c96c6', '8c6bb1',
          '88419d', '810f7c', '4d004b']
    },
    RdPu: {
      type: 'sequential',
      cbf: 42,
      3: ['fde0dd', 'fa9fb5', 'c51b8a'],
      4: ['feebe2', 'fbb4b9', 'f768a1', 'ae017e'],
      5: ['feebe2', 'fbb4b9', 'f768a1', 'c51b8a', '7a0177'],
      6: ['feebe2', 'fcc5c0', 'fa9fb5', 'f768a1', 'c51b8a', '7a0177'],
      7: ['feebe2', 'fcc5c0', 'fa9fb5', 'f768a1', 'dd3497', 'ae017e',
          '7a0177'],
      8: ['fff7f3', 'fde0dd', 'fcc5c0', 'fa9fb5', 'f768a1', 'dd3497',
          'ae017e', '7a0177'],
      9: ['fff7f3', 'fde0dd', 'fcc5c0', 'fa9fb5', 'f768a1', 'dd3497',
          'ae017e', '7a0177', '49006a']
    },
    PuRd: {
      type: 'sequential',
      cbf: 42,
      3: ['e7e1ef', 'c994c7', 'dd1c77'],
      4: ['f1eef6', 'd7b5d8', 'df65b0', 'ce1256'],
      5: ['f1eef6', 'd7b5d8', 'df65b0', 'dd1c77', '980043'],
      6: ['f1eef6', 'd4b9da', 'c994c7', 'df65b0', 'dd1c77', '980043'],
      7: ['f1eef6', 'd4b9da', 'c994c7', 'df65b0', 'e7298a', 'ce1256',
          '91003f'],
      8: ['f7f4f9', 'e7e1ef', 'd4b9da', 'c994c7', 'df65b0', 'e7298a',
          'ce1256', '91003f'],
      9: ['f7f4f9', 'e7e1ef', 'd4b9da', 'c994c7', 'df65b0', 'e7298a',
          'ce1256', '980043', '67001f']
    },
    OrRd: {
      type: 'sequential',
      cbf: 42,
      3: ['fee8c8', 'fdbb84', 'e34a33'],
      4: ['fef0d9', 'fdcc8a', 'fc8d59', 'd7301f'],
      5: ['fef0d9', 'fdcc8a', 'fc8d59', 'e34a33', 'b30000'],
      6: ['fef0d9', 'fdd49e', 'fdbb84', 'fc8d59', 'e34a33', 'b30000'],
      7: ['fef0d9', 'fdd49e', 'fdbb84', 'fc8d59', 'ef6548', 'd7301f',
          '990000'],
      8: ['fff7ec', 'fee8c8', 'fdd49e', 'fdbb84', 'fc8d59', 'ef6548',
          'd7301f', '990000'],
      9: ['fff7ec', 'fee8c8', 'fdd49e', 'fdbb84', 'fc8d59', 'ef6548',
          'd7301f', 'b30000', '7f0000']
    },
    YlOrRd: {
      type: 'sequential',
      cbf: 42,
      3: ['ffeda0', 'feb24c', 'f03b20'],
      4: ['ffffb2', 'fecc5c', 'fd8d3c', 'e31a1c'],
      5: ['ffffb2', 'fecc5c', 'fd8d3c', 'f03b20', 'bd0026'],
      6: ['ffffb2', 'fed976', 'feb24c', 'fd8d3c', 'f03b20', 'bd0026'],
      7: ['ffffb2', 'fed976', 'feb24c', 'fd8d3c', 'fc4e2a', 'e31a1c',
          'b10026'],
      8: ['ffffcc', 'ffeda0', 'fed976', 'feb24c', 'fd8d3c', 'fc4e2a',
          'e31a1c', 'b10026'],
      9: ['ffffcc', 'ffeda0', 'fed976', 'feb24c', 'fd8d3c', 'fc4e2a',
          'e31a1c', 'bd0026', '800026']
    },
    YlOrBr: {
      type: 'sequential',
      cbf: 42,
      3: ['fff7bc', 'fec44f', 'd95f0e'],
      4: ['ffffd4', 'fed98e', 'fe9929', 'cc4c02'],
      5: ['ffffd4', 'fed98e', 'fe9929', 'd95f0e', '993404'],
      6: ['ffffd4', 'fee391', 'fec44f', 'fe9929', 'd95f0e', '993404'],
      7: ['ffffd4', 'fee391', 'fec44f', 'fe9929', 'ec7014', 'cc4c02',
          '8c2d04'],
      8: ['ffffe5', 'fff7bc', 'fee391', 'fec44f', 'fe9929', 'ec7014',
          'cc4c02', '8c2d04'],
      9: ['ffffe5', 'fff7bc', 'fee391', 'fec44f', 'fe9929', 'ec7014',
          'cc4c02', '993404', '662506']
    },
    Purples: {
      type: 'sequential',
      cbf: 42,
      3: ['efedf5', 'bcbddc', '756bb1'],
      4: ['f2f0f7', 'cbc9e2', '9e9ac8', '6a51a3'],
      5: ['f2f0f7', 'cbc9e2', '9e9ac8', '756bb1', '54278f'],
      6: ['f2f0f7', 'dadaeb', 'bcbddc', '9e9ac8', '756bb1', '54278f'],
      7: ['f2f0f7', 'dadaeb', 'bcbddc', '9e9ac8', '807dba', '6a51a3',
          '4a1486'],
      8: ['fcfbfd', 'efedf5', 'dadaeb', 'bcbddc', '9e9ac8', '807dba',
          '6a51a3', '4a1486'],
      9: ['fcfbfd', 'efedf5', 'dadaeb', 'bcbddc', '9e9ac8', '807dba',
          '6a51a3', '54278f', '3f007d']
    },
    Blues: {
      type: 'sequential',
      cbf: 42,
      3: ['deebf7', '9ecae1', '3182bd'],
      4: ['eff3ff', 'bdd7e7', '6baed6', '2171b5'],
      5: ['eff3ff', 'bdd7e7', '6baed6', '3182bd', '08519c'],
      6: ['eff3ff', 'c6dbef', '9ecae1', '6baed6', '3182bd', '08519c'],
      7: ['eff3ff', 'c6dbef', '9ecae1', '6baed6', '4292c6', '2171b5',
          '084594'],
      8: ['f7fbff', 'deebf7', 'c6dbef', '9ecae1', '6baed6', '4292c6',
          '2171b5', '084594'],
      9: ['f7fbff', 'deebf7', 'c6dbef', '9ecae1', '6baed6', '4292c6',
          '2171b5', '08519c', '08306b']
    },
    Greens: {
      type: 'sequential',
      cbf: 42,
      3: ['e5f5e0', 'a1d99b', '31a354'],
      4: ['edf8e9', 'bae4b3', '74c476', '238b45'],
      5: ['edf8e9', 'bae4b3', '74c476', '31a354', '006d2c'],
      6: ['edf8e9', 'c7e9c0', 'a1d99b', '74c476', '31a354', '006d2c'],
      7: ['edf8e9', 'c7e9c0', 'a1d99b', '74c476', '41ab5d', '238b45',
          '005a32'],
      8: ['f7fcf5', 'e5f5e0', 'c7e9c0', 'a1d99b', '74c476', '41ab5d',
          '238b45', '005a32'],
      9: ['f7fcf5', 'e5f5e0', 'c7e9c0', 'a1d99b', '74c476', '41ab5d',
          '238b45', '006d2c', '00441b']
    },
    Oranges: {
      type: 'sequential',
      cbf: 42,
      3: ['fee6ce', 'fdae6b', 'e6550d'],
      4: ['feedde', 'fdbe85', 'fd8d3c', 'd94701'],
      5: ['feedde', 'fdbe85', 'fd8d3c', 'e6550d', 'a63603'],
      6: ['feedde', 'fdd0a2', 'fdae6b', 'fd8d3c', 'e6550d', 'a63603'],
      7: ['feedde', 'fdd0a2', 'fdae6b', 'fd8d3c', 'f16913', 'd94801',
          '8c2d04'],
      8: ['fff5eb', 'fee6ce', 'fdd0a2', 'fdae6b', 'fd8d3c', 'f16913',
          'd94801', '8c2d04'],
      9: ['fff5eb', 'fee6ce', 'fdd0a2', 'fdae6b', 'fd8d3c', 'f16913',
          'd94801', 'a63603', '7f2704']
    },
    Reds: {
      type: 'sequential',
      cbf: 42,
      3: ['fee0d2', 'fc9272', 'de2d26'],
      4: ['fee5d9', 'fcae91', 'fb6a4a', 'cb181d'],
      5: ['fee5d9', 'fcae91', 'fb6a4a', 'de2d26', 'a50f15'],
      6: ['fee5d9', 'fcbba1', 'fc9272', 'fb6a4a', 'de2d26', 'a50f15'],
      7: ['fee5d9', 'fcbba1', 'fc9272', 'fb6a4a', 'ef3b2c', 'cb181d',
          '99000d'],
      8: ['fff5f0', 'fee0d2', 'fcbba1', 'fc9272', 'fb6a4a', 'ef3b2c',
          'cb181d', '99000d'],
      9: ['fff5f0', 'fee0d2', 'fcbba1', 'fc9272', 'fb6a4a', 'ef3b2c',
          'cb181d', 'a50f15', '67000d']
    },
    Greys: {
      type: 'sequential',
      cbf: 42,
      3: ['f0f0f0', 'bdbdbd', '636363'],
      4: ['f7f7f7', 'cccccc', '969696', '525252'],
      5: ['f7f7f7', 'cccccc', '969696', '636363', '252525'],
      6: ['f7f7f7', 'd9d9d9', 'bdbdbd', '969696', '636363', '252525'],
      7: ['f7f7f7', 'd9d9d9', 'bdbdbd', '969696', '737373', '525252',
          '252525'],
      8: ['ffffff', 'f0f0f0', 'd9d9d9', 'bdbdbd', '969696', '737373',
          '525252', '252525'],
      9: ['ffffff', 'f0f0f0', 'd9d9d9', 'bdbdbd', '969696', '737373',
          '525252', '252525', '000000']
    },
    PuOr: {
      type: 'diverging',
      cbf: 42,
      3: ['f1a340', 'f7f7f7', '998ec3'],
      4: ['e66101', 'fdb863', 'b2abd2', '5e3c99'],
      5: ['e66101', 'fdb863', 'f7f7f7', 'b2abd2', '5e3c99'],
      6: ['b35806', 'f1a340', 'fee0b6', 'd8daeb', '998ec3', '542788'],
      7: ['b35806', 'f1a340', 'fee0b6', 'f7f7f7', 'd8daeb', '998ec3',
          '542788'],
      8: ['b35806', 'e08214', 'fdb863', 'fee0b6', 'd8daeb', 'b2abd2',
          '8073ac', '542788'],
      9: ['b35806', 'e08214', 'fdb863', 'fee0b6', 'f7f7f7', 'd8daeb',
          'b2abd2', '8073ac', '542788'],
      10: ['7f3b08', 'b35806', 'e08214', 'fdb863', 'fee0b6', 'd8daeb',
           'b2abd2', '8073ac', '542788', '2d004b'],
      11: ['7f3b08', 'b35806', 'e08214', 'fdb863', 'fee0b6', 'f7f7f7',
           'd8daeb', 'b2abd2', '8073ac', '542788', '2d004b']
    },
    BrBG: {
      type: 'diverging',
      cbf: 42,
      3: ['d8b365', 'f5f5f5', '5ab4ac'],
      4: ['a6611a', 'dfc27d', '80cdc1', '018571'],
      5: ['a6611a', 'dfc27d', 'f5f5f5', '80cdc1', '018571'],
      6: ['8c510a', 'd8b365', 'f6e8c3', 'c7eae5', '5ab4ac', '01665e'],
      7: ['8c510a', 'd8b365', 'f6e8c3', 'f5f5f5', 'c7eae5', '5ab4ac',
          '01665e'],
      8: ['8c510a', 'bf812d', 'dfc27d', 'f6e8c3', 'c7eae5', '80cdc1',
          '35978f', '01665e'],
      9: ['8c510a', 'bf812d', 'dfc27d', 'f6e8c3', 'f5f5f5', 'c7eae5',
          '80cdc1', '35978f', '01665e'],
      10: ['543005', '8c510a', 'bf812d', 'dfc27d', 'f6e8c3', 'c7eae5',
           '80cdc1', '35978f', '01665e', '003c30'],
      11: ['543005', '8c510a', 'bf812d', 'dfc27d', 'f6e8c3', 'f5f5f5',
           'c7eae5', '80cdc1', '35978f', '01665e', '003c30']
    },
    PRGn: {
      type: 'diverging',
      cbf: 42,
      3: ['af8dc3', 'f7f7f7', '7fbf7b'],
      4: ['7b3294', 'c2a5cf', 'a6dba0', '008837'],
      5: ['7b3294', 'c2a5cf', 'f7f7f7', 'a6dba0', '008837'],
      6: ['762a83', 'af8dc3', 'e7d4e8', 'd9f0d3', '7fbf7b', '1b7837'],
      7: ['762a83', 'af8dc3', 'e7d4e8', 'f7f7f7', 'd9f0d3', '7fbf7b',
          '1b7837'],
      8: ['762a83', '9970ab', 'c2a5cf', 'e7d4e8', 'd9f0d3', 'a6dba0',
          '5aae61', '1b7837'],
      9: ['762a83', '9970ab', 'c2a5cf', 'e7d4e8', 'f7f7f7', 'd9f0d3',
          'a6dba0', '5aae61', '1b7837'],
      10: ['40004b', '762a83', '9970ab', 'c2a5cf', 'e7d4e8', 'd9f0d3',
           'a6dba0', '5aae61', '1b7837', '00441b'],
      11: ['40004b', '762a83', '9970ab', 'c2a5cf', 'e7d4e8', 'f7f7f7',
           'd9f0d3', 'a6dba0', '5aae61', '1b7837', '00441b']
    },
    PiYG: {
      type: 'diverging',
      cbf: 42,
      3: ['e9a3c9', 'f7f7f7', 'a1d76a'],
      4: ['d01c8b', 'f1b6da', 'b8e186', '4dac26'],
      5: ['d01c8b', 'f1b6da', 'f7f7f7', 'b8e186', '4dac26'],
      6: ['c51b7d', 'e9a3c9', 'fde0ef', 'e6f5d0', 'a1d76a', '4d9221'],
      7: ['c51b7d', 'e9a3c9', 'fde0ef', 'f7f7f7', 'e6f5d0', 'a1d76a',
          '4d9221'],
      8: ['c51b7d', 'de77ae', 'f1b6da', 'fde0ef', 'e6f5d0', 'b8e186',
          '7fbc41', '4d9221'],
      9: ['c51b7d', 'de77ae', 'f1b6da', 'fde0ef', 'f7f7f7', 'e6f5d0',
          'b8e186', '7fbc41', '4d9221'],
      10: ['8e0152', 'c51b7d', 'de77ae', 'f1b6da', 'fde0ef', 'e6f5d0',
           'b8e186', '7fbc41', '4d9221', '276419'],
      11: ['8e0152', 'c51b7d', 'de77ae', 'f1b6da', 'fde0ef', 'f7f7f7',
           'e6f5d0', 'b8e186', '7fbc41', '4d9221', '276419']
    },
    RdBu: {
      type: 'diverging',
      cbf: 42,
      3: ['ef8a62', 'f7f7f7', '67a9cf'],
      4: ['ca0020', 'f4a582', '92c5de', '0571b0'],
      5: ['ca0020', 'f4a582', 'f7f7f7', '92c5de', '0571b0'],
      6: ['b2182b', 'ef8a62', 'fddbc7', 'd1e5f0', '67a9cf', '2166ac'],
      7: ['b2182b', 'ef8a62', 'fddbc7', 'f7f7f7', 'd1e5f0', '67a9cf',
          '2166ac'],
      8: ['b2182b', 'd6604d', 'f4a582', 'fddbc7', 'd1e5f0', '92c5de',
          '4393c3', '2166ac'],
      9: ['b2182b', 'd6604d', 'f4a582', 'fddbc7', 'f7f7f7', 'd1e5f0',
          '92c5de', '4393c3', '2166ac'],
      10: ['67001f', 'b2182b', 'd6604d', 'f4a582', 'fddbc7', 'd1e5f0',
           '92c5de', '4393c3', '2166ac', '053061'],
      11: ['67001f', 'b2182b', 'd6604d', 'f4a582', 'fddbc7', 'f7f7f7',
           'd1e5f0', '92c5de', '4393c3', '2166ac', '053061']
    },
    RdGy: {
      type: 'diverging',
      cbf: 42,
      3: ['ef8a62', 'ffffff', '999999'],
      4: ['ca0020', 'f4a582', 'bababa', '404040'],
      5: ['ca0020', 'f4a582', 'ffffff', 'bababa', '404040'],
      6: ['b2182b', 'ef8a62', 'fddbc7', 'e0e0e0', '999999', '4d4d4d'],
      7: ['b2182b', 'ef8a62', 'fddbc7', 'ffffff', 'e0e0e0', '999999',
          '4d4d4d'],
      8: ['b2182b', 'd6604d', 'f4a582', 'fddbc7', 'e0e0e0', 'bababa',
          '878787', '4d4d4d'],
      9: ['b2182b', 'd6604d', 'f4a582', 'fddbc7', 'ffffff', 'e0e0e0',
          'bababa', '878787', '4d4d4d'],
      10: ['67001f', 'b2182b', 'd6604d', 'f4a582', 'fddbc7', 'e0e0e0',
           'bababa', '878787', '4d4d4d', '1a1a1a'],
      11: ['67001f', 'b2182b', 'd6604d', 'f4a582', 'fddbc7', 'ffffff',
           'e0e0e0', 'bababa', '878787', '4d4d4d', '1a1a1a']
    },
    RdYlBu: {
      type: 'diverging',
      cbf: 42,
      3: ['fc8d59', 'ffffbf', '91bfdb'],
      4: ['d7191c', 'fdae61', 'abd9e9', '2c7bb6'],
      5: ['d7191c', 'fdae61', 'ffffbf', 'abd9e9', '2c7bb6'],
      6: ['d73027', 'fc8d59', 'fee090', 'e0f3f8', '91bfdb', '4575b4'],
      7: ['d73027', 'fc8d59', 'fee090', 'ffffbf', 'e0f3f8', '91bfdb',
          '4575b4'],
      8: ['d73027', 'f46d43', 'fdae61', 'fee090', 'e0f3f8', 'abd9e9',
          '74add1', '4575b4'],
      9: ['d73027', 'f46d43', 'fdae61', 'fee090', 'ffffbf', 'e0f3f8',
          'abd9e9', '74add1', '4575b4'],
      10: ['a50026', 'd73027', 'f46d43', 'fdae61', 'fee090', 'e0f3f8',
           'abd9e9', '74add1', '4575b4', '313695'],
      11: ['a50026', 'd73027', 'f46d43', 'fdae61', 'fee090', 'ffffbf',
           'e0f3f8', 'abd9e9', '74add1', '4575b4', '313695']
    },
    Spectral: {
      type: 'diverging',
      cbf: 0,
      3: ['fc8d59', 'ffffbf', '99d594'],
      4: ['d7191c', 'fdae61', 'abdda4', '2b83ba'],
      5: ['d7191c', 'fdae61', 'ffffbf', 'abdda4', '2b83ba'],
      6: ['d53e4f', 'fc8d59', 'fee08b', 'e6f598', '99d594', '3288bd'],
      7: ['d53e4f', 'fc8d59', 'fee08b', 'ffffbf', 'e6f598', '99d594',
          '3288bd'],
      8: ['d53e4f', 'f46d43', 'fdae61', 'fee08b', 'e6f598', 'abdda4',
          '66c2a5', '3288bd'],
      9: ['d53e4f', 'f46d43', 'fdae61', 'fee08b', 'ffffbf', 'e6f598',
          'abdda4', '66c2a5', '3288bd'],
      10: ['9e0142', 'd53e4f', 'f46d43', 'fdae61', 'fee08b', 'e6f598',
           'abdda4', '66c2a5', '3288bd', '5e4fa2'],
      11: ['9e0142', 'd53e4f', 'f46d43', 'fdae61', 'fee08b', 'ffffbf',
           'e6f598', 'abdda4', '66c2a5', '3288bd', '5e4fa2']
    },
    RdYlGn: {
      type: 'diverging',
      cbf: 0,
      3: ['fc8d59', 'ffffbf', '91cf60'],
      4: ['d7191c', 'fdae61', 'a6d96a', '1a9641'],
      5: ['d7191c', 'fdae61', 'ffffbf', 'a6d96a', '1a9641'],
      6: ['d73027', 'fc8d59', 'fee08b', 'd9ef8b', '91cf60', '1a9850'],
      7: ['d73027', 'fc8d59', 'fee08b', 'ffffbf', 'd9ef8b', '91cf60',
          '1a9850'],
      8: ['d73027', 'f46d43', 'fdae61', 'fee08b', 'd9ef8b', 'a6d96a',
          '66bd63', '1a9850'],
      9: ['d73027', 'f46d43', 'fdae61', 'fee08b', 'ffffbf', 'd9ef8b',
          'a6d96a', '66bd63', '1a9850'],
      10: ['a50026', 'd73027', 'f46d43', 'fdae61', 'fee08b', 'd9ef8b',
           'a6d96a', '66bd63', '1a9850', '006837'],
      11: ['a50026', 'd73027', 'f46d43', 'fdae61', 'fee08b', 'ffffbf',
           'd9ef8b', 'a6d96a', '66bd63', '1a9850', '006837']
    },
    Accent: {
      type: 'qualitative',
      cbf: 0,
      3: ['7fc97f', 'beaed4', 'fdc086'],
      4: ['7fc97f', 'beaed4', 'fdc086', 'ffff99'],
      5: ['7fc97f', 'beaed4', 'fdc086', 'ffff99', '386cb0'],
      6: ['7fc97f', 'beaed4', 'fdc086', 'ffff99', '386cb0', 'f0027f'],
      7: ['7fc97f', 'beaed4', 'fdc086', 'ffff99', '386cb0', 'f0027f',
          'bf5b17'],
      8: ['7fc97f', 'beaed4', 'fdc086', 'ffff99', '386cb0', 'f0027f',
          'bf5b17', '666666']
    },
    Dark2: {
      type: 'qualitative',
      cbf: 3,
      3: ['1b9e77', 'd95f02', '7570b3'],
      4: ['1b9e77', 'd95f02', '7570b3', 'e7298a'],
      5: ['1b9e77', 'd95f02', '7570b3', 'e7298a', '66a61e'],
      6: ['1b9e77', 'd95f02', '7570b3', 'e7298a', '66a61e', 'e6ab02'],
      7: ['1b9e77', 'd95f02', '7570b3', 'e7298a', '66a61e', 'e6ab02',
          'a6761d'],
      8: ['1b9e77', 'd95f02', '7570b3', 'e7298a', '66a61e', 'e6ab02',
          'a6761d', '666666']
    },
    Paired: {
      type: 'qualitative',
      cbf: 4,
      3: ['a6cee3', '1f78b4', 'b2df8a'],
      4: ['a6cee3', '1f78b4', 'b2df8a', '33a02c'],
      5: ['a6cee3', '1f78b4', 'b2df8a', '33a02c', 'fb9a99'],
      6: ['a6cee3', '1f78b4', 'b2df8a', '33a02c', 'fb9a99', 'e31a1c'],
      7: ['a6cee3', '1f78b4', 'b2df8a', '33a02c', 'fb9a99', 'e31a1c',
          'fdbf6f'],
      8: ['a6cee3', '1f78b4', 'b2df8a', '33a02c', 'fb9a99', 'e31a1c',
          'fdbf6f', 'ff7f00'],
      9: ['a6cee3', '1f78b4', 'b2df8a', '33a02c', 'fb9a99', 'e31a1c',
          'fdbf6f', 'ff7f00', 'cab2d6'],
      10: ['a6cee3', '1f78b4', 'b2df8a', '33a02c', 'fb9a99', 'e31a1c',
           'fdbf6f', 'ff7f00', 'cab2d6', '6a3d9a'],
      11: ['a6cee3', '1f78b4', 'b2df8a', '33a02c', 'fb9a99', 'e31a1c',
           'fdbf6f', 'ff7f00', 'cab2d6', '6a3d9a', 'ffff99'],
      12: ['a6cee3', '1f78b4', 'b2df8a', '33a02c', 'fb9a99', 'e31a1c',
           'fdbf6f', 'ff7f00', 'cab2d6', '6a3d9a', 'ffff99', 'b15928']
    },
    Pastel1: {
      type: 'qualitative',
      cbf: 0,
      3: ['fbb4ae', 'b3cde3', 'ccebc5'],
      4: ['fbb4ae', 'b3cde3', 'ccebc5', 'decbe4'],
      5: ['fbb4ae', 'b3cde3', 'ccebc5', 'decbe4', 'fed9a6'],
      6: ['fbb4ae', 'b3cde3', 'ccebc5', 'decbe4', 'fed9a6', 'ffffcc'],
      7: ['fbb4ae', 'b3cde3', 'ccebc5', 'decbe4', 'fed9a6', 'ffffcc',
          'e5d8bd'],
      8: ['fbb4ae', 'b3cde3', 'ccebc5', 'decbe4', 'fed9a6', 'ffffcc',
          'e5d8bd', 'fddaec'],
      9: ['fbb4ae', 'b3cde3', 'ccebc5', 'decbe4', 'fed9a6', 'ffffcc',
          'e5d8bd', 'fddaec', 'f2f2f2']
    },
    Pastel2: {
      type: 'qualitative',
      cbf: 0,
      3: ['b3e2cd', 'fdcdac', 'cbd5e8'],
      4: ['b3e2cd', 'fdcdac', 'cbd5e8', 'f4cae4'],
      5: ['b3e2cd', 'fdcdac', 'cbd5e8', 'f4cae4', 'e6f5c9'],
      6: ['b3e2cd', 'fdcdac', 'cbd5e8', 'f4cae4', 'e6f5c9', 'fff2ae'],
      7: ['b3e2cd', 'fdcdac', 'cbd5e8', 'f4cae4', 'e6f5c9', 'fff2ae',
          'f1e2cc'],
      8: ['b3e2cd', 'fdcdac', 'cbd5e8', 'f4cae4', 'e6f5c9', 'fff2ae',
          'f1e2cc', 'cccccc']
    },
    Set1: {
      type: 'qualitative',
      cbf: 0,
      3: ['e41a1c', '377eb8', '4daf4a'],
      4: ['e41a1c', '377eb8', '4daf4a', '984ea3'],
      5: ['e41a1c', '377eb8', '4daf4a', '984ea3', 'ff7f00'],
      6: ['e41a1c', '377eb8', '4daf4a', '984ea3', 'ff7f00', 'ffff33'],
      7: ['e41a1c', '377eb8', '4daf4a', '984ea3', 'ff7f00', 'ffff33',
          'a65628'],
      8: ['e41a1c', '377eb8', '4daf4a', '984ea3', 'ff7f00', 'ffff33',
          'a65628', 'f781bf'],
      9: ['e41a1c', '377eb8', '4daf4a', '984ea3', 'ff7f00', 'ffff33',
          'a65628', 'f781bf', '999999']
    },
    Set2: {
      type: 'qualitative',
      cbf: 3,
      3: ['66c2a5', 'fc8d62', '8da0cb'],
      4: ['66c2a5', 'fc8d62', '8da0cb', 'e78ac3'],
      5: ['66c2a5', 'fc8d62', '8da0cb', 'e78ac3', 'a6d854'],
      6: ['66c2a5', 'fc8d62', '8da0cb', 'e78ac3', 'a6d854', 'ffd92f'],
      7: ['66c2a5', 'fc8d62', '8da0cb', 'e78ac3', 'a6d854', 'ffd92f',
          'e5c494'],
      8: ['66c2a5', 'fc8d62', '8da0cb', 'e78ac3', 'a6d854', 'ffd92f',
          'e5c494', 'b3b3b3']
    },
    Set3: {
      type: 'qualitative',
      cbf: 0,
      3: ['8dd3c7', 'ffffb3', 'bebada'],
      4: ['8dd3c7', 'ffffb3', 'bebada', 'fb8072'],
      5: ['8dd3c7', 'ffffb3', 'bebada', 'fb8072', '80b1d3'],
      6: ['8dd3c7', 'ffffb3', 'bebada', 'fb8072', '80b1d3', 'fdb462'],
      7: ['8dd3c7', 'ffffb3', 'bebada', 'fb8072', '80b1d3', 'fdb462',
          'b3de69'],
      8: ['8dd3c7', 'ffffb3', 'bebada', 'fb8072', '80b1d3', 'fdb462',
          'b3de69', 'fccde5'],
      9: ['8dd3c7', 'ffffb3', 'bebada', 'fb8072', '80b1d3', 'fdb462',
          'b3de69', 'fccde5', 'd9d9d9'],
      10: ['8dd3c7', 'ffffb3', 'bebada', 'fb8072', '80b1d3', 'fdb462',
           'b3de69', 'fccde5', 'd9d9d9', 'bc80bd'],
      11: ['8dd3c7', 'ffffb3', 'bebada', 'fb8072', '80b1d3', 'fdb462',
           'b3de69', 'fccde5', 'd9d9d9', 'bc80bd', 'ccebc5'],
      12: ['8dd3c7', 'ffffb3', 'bebada', 'fb8072', '80b1d3', 'fdb462',
           'b3de69', 'fccde5', 'd9d9d9', 'bc80bd', 'ccebc5', 'ffed6f']
    }
  };

  for (var name in schemes) {
    var scheme = schemes[name];
    scheme = palette.Scheme.fromPalettes(
      'cb-' + name, [scheme.type, 'cb-' + scheme.type], scheme, 12, scheme.cbf);
    palette.register(scheme);
  }
})();


! function() {
 "use strict";
 "SVGPathSeg" in window || (window.SVGPathSeg = function(t, e, i) {
  this.pathSegType = t, this.pathSegTypeAsLetter = e, this._owningPathSegList = i
 }, window.SVGPathSeg.prototype.classname = "SVGPathSeg", window.SVGPathSeg.PATHSEG_UNKNOWN = 0, window.SVGPathSeg.PATHSEG_CLOSEPATH = 1, window.SVGPathSeg.PATHSEG_MOVETO_ABS = 2, window.SVGPathSeg.PATHSEG_MOVETO_REL = 3, window.SVGPathSeg.PATHSEG_LINETO_ABS = 4, window.SVGPathSeg.PATHSEG_LINETO_REL = 5, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_ABS = 6, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_REL = 7, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_ABS = 8, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_REL = 9, window.SVGPathSeg.PATHSEG_ARC_ABS = 10, window.SVGPathSeg.PATHSEG_ARC_REL = 11, window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_ABS = 12, window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_REL = 13, window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_ABS = 14, window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_REL = 15, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_ABS = 16, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_REL = 17, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_ABS = 18, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_REL = 19, window.SVGPathSeg.prototype._segmentChanged = function() {
  this._owningPathSegList && this._owningPathSegList.segmentChanged(this)
 }, window.SVGPathSegClosePath = function(t) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CLOSEPATH, "z", t)
 }, window.SVGPathSegClosePath.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegClosePath.prototype.toString = function() {
  return "[object SVGPathSegClosePath]"
 }, window.SVGPathSegClosePath.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter
 }, window.SVGPathSegClosePath.prototype.clone = function() {
  return new window.SVGPathSegClosePath(void 0)
 }, window.SVGPathSegMovetoAbs = function(t, e, i) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_MOVETO_ABS, "M", t), this._x = e, this._y = i
 }, window.SVGPathSegMovetoAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegMovetoAbs.prototype.toString = function() {
  return "[object SVGPathSegMovetoAbs]"
 }, window.SVGPathSegMovetoAbs.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x + " " + this._y
 }, window.SVGPathSegMovetoAbs.prototype.clone = function() {
  return new window.SVGPathSegMovetoAbs(void 0, this._x, this._y)
 }, Object.defineProperty(window.SVGPathSegMovetoAbs.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegMovetoAbs.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegMovetoRel = function(t, e, i) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_MOVETO_REL, "m", t), this._x = e, this._y = i
 }, window.SVGPathSegMovetoRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegMovetoRel.prototype.toString = function() {
  return "[object SVGPathSegMovetoRel]"
 }, window.SVGPathSegMovetoRel.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x + " " + this._y
 }, window.SVGPathSegMovetoRel.prototype.clone = function() {
  return new window.SVGPathSegMovetoRel(void 0, this._x, this._y)
 }, Object.defineProperty(window.SVGPathSegMovetoRel.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegMovetoRel.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegLinetoAbs = function(t, e, i) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_LINETO_ABS, "L", t), this._x = e, this._y = i
 }, window.SVGPathSegLinetoAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegLinetoAbs.prototype.toString = function() {
  return "[object SVGPathSegLinetoAbs]"
 }, window.SVGPathSegLinetoAbs.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x + " " + this._y
 }, window.SVGPathSegLinetoAbs.prototype.clone = function() {
  return new window.SVGPathSegLinetoAbs(void 0, this._x, this._y)
 }, Object.defineProperty(window.SVGPathSegLinetoAbs.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegLinetoAbs.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegLinetoRel = function(t, e, i) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_LINETO_REL, "l", t), this._x = e, this._y = i
 }, window.SVGPathSegLinetoRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegLinetoRel.prototype.toString = function() {
  return "[object SVGPathSegLinetoRel]"
 }, window.SVGPathSegLinetoRel.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x + " " + this._y
 }, window.SVGPathSegLinetoRel.prototype.clone = function() {
  return new window.SVGPathSegLinetoRel(void 0, this._x, this._y)
 }, Object.defineProperty(window.SVGPathSegLinetoRel.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegLinetoRel.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegCurvetoCubicAbs = function(t, e, i, n, r, o, s) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_ABS, "C", t), this._x = e, this._y = i, this._x1 = n, this._y1 = r, this._x2 = o, this._y2 = s
 }, window.SVGPathSegCurvetoCubicAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoCubicAbs.prototype.toString = function() {
  return "[object SVGPathSegCurvetoCubicAbs]"
 }, window.SVGPathSegCurvetoCubicAbs.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x1 + " " + this._y1 + " " + this._x2 + " " + this._y2 + " " + this._x + " " + this._y
 }, window.SVGPathSegCurvetoCubicAbs.prototype.clone = function() {
  return new window.SVGPathSegCurvetoCubicAbs(void 0, this._x, this._y, this._x1, this._y1, this._x2, this._y2)
 }, Object.defineProperty(window.SVGPathSegCurvetoCubicAbs.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicAbs.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicAbs.prototype, "x1", {
  get: function() {
   return this._x1
  },
  set: function(t) {
   this._x1 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicAbs.prototype, "y1", {
  get: function() {
   return this._y1
  },
  set: function(t) {
   this._y1 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicAbs.prototype, "x2", {
  get: function() {
   return this._x2
  },
  set: function(t) {
   this._x2 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicAbs.prototype, "y2", {
  get: function() {
   return this._y2
  },
  set: function(t) {
   this._y2 = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegCurvetoCubicRel = function(t, e, i, n, r, o, s) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_REL, "c", t), this._x = e, this._y = i, this._x1 = n, this._y1 = r, this._x2 = o, this._y2 = s
 }, window.SVGPathSegCurvetoCubicRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoCubicRel.prototype.toString = function() {
  return "[object SVGPathSegCurvetoCubicRel]"
 }, window.SVGPathSegCurvetoCubicRel.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x1 + " " + this._y1 + " " + this._x2 + " " + this._y2 + " " + this._x + " " + this._y
 }, window.SVGPathSegCurvetoCubicRel.prototype.clone = function() {
  return new window.SVGPathSegCurvetoCubicRel(void 0, this._x, this._y, this._x1, this._y1, this._x2, this._y2)
 }, Object.defineProperty(window.SVGPathSegCurvetoCubicRel.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicRel.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicRel.prototype, "x1", {
  get: function() {
   return this._x1
  },
  set: function(t) {
   this._x1 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicRel.prototype, "y1", {
  get: function() {
   return this._y1
  },
  set: function(t) {
   this._y1 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicRel.prototype, "x2", {
  get: function() {
   return this._x2
  },
  set: function(t) {
   this._x2 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicRel.prototype, "y2", {
  get: function() {
   return this._y2
  },
  set: function(t) {
   this._y2 = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegCurvetoQuadraticAbs = function(t, e, i, n, r) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_ABS, "Q", t), this._x = e, this._y = i, this._x1 = n, this._y1 = r
 }, window.SVGPathSegCurvetoQuadraticAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoQuadraticAbs.prototype.toString = function() {
  return "[object SVGPathSegCurvetoQuadraticAbs]"
 }, window.SVGPathSegCurvetoQuadraticAbs.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x1 + " " + this._y1 + " " + this._x + " " + this._y
 }, window.SVGPathSegCurvetoQuadraticAbs.prototype.clone = function() {
  return new window.SVGPathSegCurvetoQuadraticAbs(void 0, this._x, this._y, this._x1, this._y1)
 }, Object.defineProperty(window.SVGPathSegCurvetoQuadraticAbs.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticAbs.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticAbs.prototype, "x1", {
  get: function() {
   return this._x1
  },
  set: function(t) {
   this._x1 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticAbs.prototype, "y1", {
  get: function() {
   return this._y1
  },
  set: function(t) {
   this._y1 = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegCurvetoQuadraticRel = function(t, e, i, n, r) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_REL, "q", t), this._x = e, this._y = i, this._x1 = n, this._y1 = r
 }, window.SVGPathSegCurvetoQuadraticRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoQuadraticRel.prototype.toString = function() {
  return "[object SVGPathSegCurvetoQuadraticRel]"
 }, window.SVGPathSegCurvetoQuadraticRel.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x1 + " " + this._y1 + " " + this._x + " " + this._y
 }, window.SVGPathSegCurvetoQuadraticRel.prototype.clone = function() {
  return new window.SVGPathSegCurvetoQuadraticRel(void 0, this._x, this._y, this._x1, this._y1)
 }, Object.defineProperty(window.SVGPathSegCurvetoQuadraticRel.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticRel.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticRel.prototype, "x1", {
  get: function() {
   return this._x1
  },
  set: function(t) {
   this._x1 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticRel.prototype, "y1", {
  get: function() {
   return this._y1
  },
  set: function(t) {
   this._y1 = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegArcAbs = function(t, e, i, n, r, o, s, a) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_ARC_ABS, "A", t), this._x = e, this._y = i, this._r1 = n, this._r2 = r, this._angle = o, this._largeArcFlag = s, this._sweepFlag = a
 }, window.SVGPathSegArcAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegArcAbs.prototype.toString = function() {
  return "[object SVGPathSegArcAbs]"
 }, window.SVGPathSegArcAbs.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._r1 + " " + this._r2 + " " + this._angle + " " + (this._largeArcFlag ? "1" : "0") + " " + (this._sweepFlag ? "1" : "0") + " " + this._x + " " + this._y
 }, window.SVGPathSegArcAbs.prototype.clone = function() {
  return new window.SVGPathSegArcAbs(void 0, this._x, this._y, this._r1, this._r2, this._angle, this._largeArcFlag, this._sweepFlag)
 }, Object.defineProperty(window.SVGPathSegArcAbs.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegArcAbs.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegArcAbs.prototype, "r1", {
  get: function() {
   return this._r1
  },
  set: function(t) {
   this._r1 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegArcAbs.prototype, "r2", {
  get: function() {
   return this._r2
  },
  set: function(t) {
   this._r2 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegArcAbs.prototype, "angle", {
  get: function() {
   return this._angle
  },
  set: function(t) {
   this._angle = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegArcAbs.prototype, "largeArcFlag", {
  get: function() {
   return this._largeArcFlag
  },
  set: function(t) {
   this._largeArcFlag = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegArcAbs.prototype, "sweepFlag", {
  get: function() {
   return this._sweepFlag
  },
  set: function(t) {
   this._sweepFlag = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegArcRel = function(t, e, i, n, r, o, s, a) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_ARC_REL, "a", t), this._x = e, this._y = i, this._r1 = n, this._r2 = r, this._angle = o, this._largeArcFlag = s, this._sweepFlag = a
 }, window.SVGPathSegArcRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegArcRel.prototype.toString = function() {
  return "[object SVGPathSegArcRel]"
 }, window.SVGPathSegArcRel.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._r1 + " " + this._r2 + " " + this._angle + " " + (this._largeArcFlag ? "1" : "0") + " " + (this._sweepFlag ? "1" : "0") + " " + this._x + " " + this._y
 }, window.SVGPathSegArcRel.prototype.clone = function() {
  return new window.SVGPathSegArcRel(void 0, this._x, this._y, this._r1, this._r2, this._angle, this._largeArcFlag, this._sweepFlag)
 }, Object.defineProperty(window.SVGPathSegArcRel.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegArcRel.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegArcRel.prototype, "r1", {
  get: function() {
   return this._r1
  },
  set: function(t) {
   this._r1 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegArcRel.prototype, "r2", {
  get: function() {
   return this._r2
  },
  set: function(t) {
   this._r2 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegArcRel.prototype, "angle", {
  get: function() {
   return this._angle
  },
  set: function(t) {
   this._angle = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegArcRel.prototype, "largeArcFlag", {
  get: function() {
   return this._largeArcFlag
  },
  set: function(t) {
   this._largeArcFlag = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegArcRel.prototype, "sweepFlag", {
  get: function() {
   return this._sweepFlag
  },
  set: function(t) {
   this._sweepFlag = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegLinetoHorizontalAbs = function(t, e) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_ABS, "H", t), this._x = e
 }, window.SVGPathSegLinetoHorizontalAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegLinetoHorizontalAbs.prototype.toString = function() {
  return "[object SVGPathSegLinetoHorizontalAbs]"
 }, window.SVGPathSegLinetoHorizontalAbs.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x
 }, window.SVGPathSegLinetoHorizontalAbs.prototype.clone = function() {
  return new window.SVGPathSegLinetoHorizontalAbs(void 0, this._x)
 }, Object.defineProperty(window.SVGPathSegLinetoHorizontalAbs.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegLinetoHorizontalRel = function(t, e) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_REL, "h", t), this._x = e
 }, window.SVGPathSegLinetoHorizontalRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegLinetoHorizontalRel.prototype.toString = function() {
  return "[object SVGPathSegLinetoHorizontalRel]"
 }, window.SVGPathSegLinetoHorizontalRel.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x
 }, window.SVGPathSegLinetoHorizontalRel.prototype.clone = function() {
  return new window.SVGPathSegLinetoHorizontalRel(void 0, this._x)
 }, Object.defineProperty(window.SVGPathSegLinetoHorizontalRel.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegLinetoVerticalAbs = function(t, e) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_ABS, "V", t), this._y = e
 }, window.SVGPathSegLinetoVerticalAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegLinetoVerticalAbs.prototype.toString = function() {
  return "[object SVGPathSegLinetoVerticalAbs]"
 }, window.SVGPathSegLinetoVerticalAbs.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._y
 }, window.SVGPathSegLinetoVerticalAbs.prototype.clone = function() {
  return new window.SVGPathSegLinetoVerticalAbs(void 0, this._y)
 }, Object.defineProperty(window.SVGPathSegLinetoVerticalAbs.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegLinetoVerticalRel = function(t, e) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_REL, "v", t), this._y = e
 }, window.SVGPathSegLinetoVerticalRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegLinetoVerticalRel.prototype.toString = function() {
  return "[object SVGPathSegLinetoVerticalRel]"
 }, window.SVGPathSegLinetoVerticalRel.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._y
 }, window.SVGPathSegLinetoVerticalRel.prototype.clone = function() {
  return new window.SVGPathSegLinetoVerticalRel(void 0, this._y)
 }, Object.defineProperty(window.SVGPathSegLinetoVerticalRel.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegCurvetoCubicSmoothAbs = function(t, e, i, n, r) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_ABS, "S", t), this._x = e, this._y = i, this._x2 = n, this._y2 = r
 }, window.SVGPathSegCurvetoCubicSmoothAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoCubicSmoothAbs.prototype.toString = function() {
  return "[object SVGPathSegCurvetoCubicSmoothAbs]"
 }, window.SVGPathSegCurvetoCubicSmoothAbs.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x2 + " " + this._y2 + " " + this._x + " " + this._y
 }, window.SVGPathSegCurvetoCubicSmoothAbs.prototype.clone = function() {
  return new window.SVGPathSegCurvetoCubicSmoothAbs(void 0, this._x, this._y, this._x2, this._y2)
 }, Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothAbs.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothAbs.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothAbs.prototype, "x2", {
  get: function() {
   return this._x2
  },
  set: function(t) {
   this._x2 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothAbs.prototype, "y2", {
  get: function() {
   return this._y2
  },
  set: function(t) {
   this._y2 = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegCurvetoCubicSmoothRel = function(t, e, i, n, r) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_REL, "s", t), this._x = e, this._y = i, this._x2 = n, this._y2 = r
 }, window.SVGPathSegCurvetoCubicSmoothRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoCubicSmoothRel.prototype.toString = function() {
  return "[object SVGPathSegCurvetoCubicSmoothRel]"
 }, window.SVGPathSegCurvetoCubicSmoothRel.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x2 + " " + this._y2 + " " + this._x + " " + this._y
 }, window.SVGPathSegCurvetoCubicSmoothRel.prototype.clone = function() {
  return new window.SVGPathSegCurvetoCubicSmoothRel(void 0, this._x, this._y, this._x2, this._y2)
 }, Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothRel.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothRel.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothRel.prototype, "x2", {
  get: function() {
   return this._x2
  },
  set: function(t) {
   this._x2 = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoCubicSmoothRel.prototype, "y2", {
  get: function() {
   return this._y2
  },
  set: function(t) {
   this._y2 = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegCurvetoQuadraticSmoothAbs = function(t, e, i) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_ABS, "T", t), this._x = e, this._y = i
 }, window.SVGPathSegCurvetoQuadraticSmoothAbs.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoQuadraticSmoothAbs.prototype.toString = function() {
  return "[object SVGPathSegCurvetoQuadraticSmoothAbs]"
 }, window.SVGPathSegCurvetoQuadraticSmoothAbs.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x + " " + this._y
 }, window.SVGPathSegCurvetoQuadraticSmoothAbs.prototype.clone = function() {
  return new window.SVGPathSegCurvetoQuadraticSmoothAbs(void 0, this._x, this._y)
 }, Object.defineProperty(window.SVGPathSegCurvetoQuadraticSmoothAbs.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticSmoothAbs.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathSegCurvetoQuadraticSmoothRel = function(t, e, i) {
  window.SVGPathSeg.call(this, window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_REL, "t", t), this._x = e, this._y = i
 }, window.SVGPathSegCurvetoQuadraticSmoothRel.prototype = Object.create(window.SVGPathSeg.prototype), window.SVGPathSegCurvetoQuadraticSmoothRel.prototype.toString = function() {
  return "[object SVGPathSegCurvetoQuadraticSmoothRel]"
 }, window.SVGPathSegCurvetoQuadraticSmoothRel.prototype._asPathString = function() {
  return this.pathSegTypeAsLetter + " " + this._x + " " + this._y
 }, window.SVGPathSegCurvetoQuadraticSmoothRel.prototype.clone = function() {
  return new window.SVGPathSegCurvetoQuadraticSmoothRel(void 0, this._x, this._y)
 }, Object.defineProperty(window.SVGPathSegCurvetoQuadraticSmoothRel.prototype, "x", {
  get: function() {
   return this._x
  },
  set: function(t) {
   this._x = t, this._segmentChanged()
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathSegCurvetoQuadraticSmoothRel.prototype, "y", {
  get: function() {
   return this._y
  },
  set: function(t) {
   this._y = t, this._segmentChanged()
  },
  enumerable: !0
 }), window.SVGPathElement.prototype.createSVGPathSegClosePath = function() {
  return new window.SVGPathSegClosePath(void 0)
 }, window.SVGPathElement.prototype.createSVGPathSegMovetoAbs = function(t, e) {
  return new window.SVGPathSegMovetoAbs(void 0, t, e)
 }, window.SVGPathElement.prototype.createSVGPathSegMovetoRel = function(t, e) {
  return new window.SVGPathSegMovetoRel(void 0, t, e)
 }, window.SVGPathElement.prototype.createSVGPathSegLinetoAbs = function(t, e) {
  return new window.SVGPathSegLinetoAbs(void 0, t, e)
 }, window.SVGPathElement.prototype.createSVGPathSegLinetoRel = function(t, e) {
  return new window.SVGPathSegLinetoRel(void 0, t, e)
 }, window.SVGPathElement.prototype.createSVGPathSegCurvetoCubicAbs = function(t, e, i, n, r, o) {
  return new window.SVGPathSegCurvetoCubicAbs(void 0, t, e, i, n, r, o)
 }, window.SVGPathElement.prototype.createSVGPathSegCurvetoCubicRel = function(t, e, i, n, r, o) {
  return new window.SVGPathSegCurvetoCubicRel(void 0, t, e, i, n, r, o)
 }, window.SVGPathElement.prototype.createSVGPathSegCurvetoQuadraticAbs = function(t, e, i, n) {
  return new window.SVGPathSegCurvetoQuadraticAbs(void 0, t, e, i, n)
 }, window.SVGPathElement.prototype.createSVGPathSegCurvetoQuadraticRel = function(t, e, i, n) {
  return new window.SVGPathSegCurvetoQuadraticRel(void 0, t, e, i, n)
 }, window.SVGPathElement.prototype.createSVGPathSegArcAbs = function(t, e, i, n, r, o, s) {
  return new window.SVGPathSegArcAbs(void 0, t, e, i, n, r, o, s)
 }, window.SVGPathElement.prototype.createSVGPathSegArcRel = function(t, e, i, n, r, o, s) {
  return new window.SVGPathSegArcRel(void 0, t, e, i, n, r, o, s)
 }, window.SVGPathElement.prototype.createSVGPathSegLinetoHorizontalAbs = function(t) {
  return new window.SVGPathSegLinetoHorizontalAbs(void 0, t)
 }, window.SVGPathElement.prototype.createSVGPathSegLinetoHorizontalRel = function(t) {
  return new window.SVGPathSegLinetoHorizontalRel(void 0, t)
 }, window.SVGPathElement.prototype.createSVGPathSegLinetoVerticalAbs = function(t) {
  return new window.SVGPathSegLinetoVerticalAbs(void 0, t)
 }, window.SVGPathElement.prototype.createSVGPathSegLinetoVerticalRel = function(t) {
  return new window.SVGPathSegLinetoVerticalRel(void 0, t)
 }, window.SVGPathElement.prototype.createSVGPathSegCurvetoCubicSmoothAbs = function(t, e, i, n) {
  return new window.SVGPathSegCurvetoCubicSmoothAbs(void 0, t, e, i, n)
 }, window.SVGPathElement.prototype.createSVGPathSegCurvetoCubicSmoothRel = function(t, e, i, n) {
  return new window.SVGPathSegCurvetoCubicSmoothRel(void 0, t, e, i, n)
 }, window.SVGPathElement.prototype.createSVGPathSegCurvetoQuadraticSmoothAbs = function(t, e) {
  return new window.SVGPathSegCurvetoQuadraticSmoothAbs(void 0, t, e)
 }, window.SVGPathElement.prototype.createSVGPathSegCurvetoQuadraticSmoothRel = function(t, e) {
  return new window.SVGPathSegCurvetoQuadraticSmoothRel(void 0, t, e)
 }, "getPathSegAtLength" in window.SVGPathElement.prototype || (window.SVGPathElement.prototype.getPathSegAtLength = function(t) {
  if (void 0 === t || !isFinite(t)) throw "Invalid arguments.";
  var e = document.createElementNS("http://www.w3.org/2000/svg", "path");
  e.setAttribute("d", this.getAttribute("d"));
  var i = e.pathSegList.numberOfItems - 1;
  if (i <= 0) return 0;
  do {
   if (e.pathSegList.removeItem(i), t > e.getTotalLength()) break;
   i--
  } while (0 < i);
  return i
 })), "SVGPathSegList" in window && "appendItem" in window.SVGPathSegList.prototype || (window.SVGPathSegList = function(t) {
  this._pathElement = t, this._list = this._parsePath(this._pathElement.getAttribute("d")), this._mutationObserverConfig = {
   attributes: !0,
   attributeFilter: ["d"]
  }, this._pathElementMutationObserver = new MutationObserver(this._updateListFromPathMutations.bind(this)), this._pathElementMutationObserver.observe(this._pathElement, this._mutationObserverConfig)
 }, window.SVGPathSegList.prototype.classname = "SVGPathSegList", Object.defineProperty(window.SVGPathSegList.prototype, "numberOfItems", {
  get: function() {
   return this._checkPathSynchronizedToList(), this._list.length
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathElement.prototype, "pathSegList", {
  get: function() {
   return this._pathSegList || (this._pathSegList = new window.SVGPathSegList(this)), this._pathSegList
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathElement.prototype, "normalizedPathSegList", {
  get: function() {
   return this.pathSegList
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathElement.prototype, "animatedPathSegList", {
  get: function() {
   return this.pathSegList
  },
  enumerable: !0
 }), Object.defineProperty(window.SVGPathElement.prototype, "animatedNormalizedPathSegList", {
  get: function() {
   return this.pathSegList
  },
  enumerable: !0
 }), window.SVGPathSegList.prototype._checkPathSynchronizedToList = function() {
  this._updateListFromPathMutations(this._pathElementMutationObserver.takeRecords())
 }, window.SVGPathSegList.prototype._updateListFromPathMutations = function(t) {
  if (this._pathElement) {
   var e = !1;
   t.forEach(function(t) {
    "d" == t.attributeName && (e = !0)
   }), e && (this._list = this._parsePath(this._pathElement.getAttribute("d")))
  }
 }, window.SVGPathSegList.prototype._writeListToPath = function() {
  this._pathElementMutationObserver.disconnect(), this._pathElement.setAttribute("d", window.SVGPathSegList._pathSegArrayAsString(this._list)), this._pathElementMutationObserver.observe(this._pathElement, this._mutationObserverConfig)
 }, window.SVGPathSegList.prototype.segmentChanged = function(t) {
  this._writeListToPath()
 }, window.SVGPathSegList.prototype.clear = function() {
  this._checkPathSynchronizedToList(), this._list.forEach(function(t) {
   t._owningPathSegList = null
  }), this._list = [], this._writeListToPath()
 }, window.SVGPathSegList.prototype.initialize = function(t) {
  return this._checkPathSynchronizedToList(), this._list = [t], (t._owningPathSegList = this)._writeListToPath(), t
 }, window.SVGPathSegList.prototype._checkValidIndex = function(t) {
  if (isNaN(t) || t < 0 || t >= this.numberOfItems) throw "INDEX_SIZE_ERR"
 }, window.SVGPathSegList.prototype.getItem = function(t) {
  return this._checkPathSynchronizedToList(), this._checkValidIndex(t), this._list[t]
 }, window.SVGPathSegList.prototype.insertItemBefore = function(t, e) {
  return this._checkPathSynchronizedToList(), e > this.numberOfItems && (e = this.numberOfItems), t._owningPathSegList && (t = t.clone()), this._list.splice(e, 0, t), (t._owningPathSegList = this)._writeListToPath(), t
 }, window.SVGPathSegList.prototype.replaceItem = function(t, e) {
  return this._checkPathSynchronizedToList(), t._owningPathSegList && (t = t.clone()), this._checkValidIndex(e), ((this._list[e] = t)._owningPathSegList = this)._writeListToPath(), t
 }, window.SVGPathSegList.prototype.removeItem = function(t) {
  this._checkPathSynchronizedToList(), this._checkValidIndex(t);
  var e = this._list[t];
  return this._list.splice(t, 1), this._writeListToPath(), e
 }, window.SVGPathSegList.prototype.appendItem = function(t) {
  return this._checkPathSynchronizedToList(), t._owningPathSegList && (t = t.clone()), this._list.push(t), (t._owningPathSegList = this)._writeListToPath(), t
 }, window.SVGPathSegList._pathSegArrayAsString = function(t) {
  var e = "",
   i = !0;
  return t.forEach(function(t) {
   i ? (i = !1, e += t._asPathString()) : e += " " + t._asPathString()
  }), e
 }, window.SVGPathSegList.prototype._parsePath = function(t) {
  if (!t || 0 == t.length) return [];

  function e() {
   this.pathSegList = []
  }
  var n = this;

  function i(t) {
   this._string = t, this._currentIndex = 0, this._endIndex = this._string.length, this._previousCommand = window.SVGPathSeg.PATHSEG_UNKNOWN, this._skipOptionalSpaces()
  }
  e.prototype.appendSegment = function(t) {
   this.pathSegList.push(t)
  }, i.prototype._isCurrentSpace = function() {
   var t = this._string[this._currentIndex];
   return t <= " " && (" " == t || "\n" == t || "\t" == t || "\r" == t || "\f" == t)
  }, i.prototype._skipOptionalSpaces = function() {
   for (; this._currentIndex < this._endIndex && this._isCurrentSpace();) this._currentIndex++;
   return this._currentIndex < this._endIndex
  }, i.prototype._skipOptionalSpacesOrDelimiter = function() {
   return !(this._currentIndex < this._endIndex && !this._isCurrentSpace() && "," != this._string.charAt(this._currentIndex)) && (this._skipOptionalSpaces() && this._currentIndex < this._endIndex && "," == this._string.charAt(this._currentIndex) && (this._currentIndex++, this._skipOptionalSpaces()), this._currentIndex < this._endIndex)
  }, i.prototype.hasMoreData = function() {
   return this._currentIndex < this._endIndex
  }, i.prototype.peekSegmentType = function() {
   var t = this._string[this._currentIndex];
   return this._pathSegTypeFromChar(t)
  }, i.prototype._pathSegTypeFromChar = function(t) {
   switch (t) {
    case "Z":
    case "z":
     return window.SVGPathSeg.PATHSEG_CLOSEPATH;
    case "M":
     return window.SVGPathSeg.PATHSEG_MOVETO_ABS;
    case "m":
     return window.SVGPathSeg.PATHSEG_MOVETO_REL;
    case "L":
     return window.SVGPathSeg.PATHSEG_LINETO_ABS;
    case "l":
     return window.SVGPathSeg.PATHSEG_LINETO_REL;
    case "C":
     return window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_ABS;
    case "c":
     return window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_REL;
    case "Q":
     return window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_ABS;
    case "q":
     return window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_REL;
    case "A":
     return window.SVGPathSeg.PATHSEG_ARC_ABS;
    case "a":
     return window.SVGPathSeg.PATHSEG_ARC_REL;
    case "H":
     return window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_ABS;
    case "h":
     return window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_REL;
    case "V":
     return window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_ABS;
    case "v":
     return window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_REL;
    case "S":
     return window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_ABS;
    case "s":
     return window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_REL;
    case "T":
     return window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_ABS;
    case "t":
     return window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_REL;
    default:
     return window.SVGPathSeg.PATHSEG_UNKNOWN
   }
  }, i.prototype._nextCommandHelper = function(t, e) {
   return ("+" == t || "-" == t || "." == t || "0" <= t && t <= "9") && e != window.SVGPathSeg.PATHSEG_CLOSEPATH ? e == window.SVGPathSeg.PATHSEG_MOVETO_ABS ? window.SVGPathSeg.PATHSEG_LINETO_ABS : e == window.SVGPathSeg.PATHSEG_MOVETO_REL ? window.SVGPathSeg.PATHSEG_LINETO_REL : e : window.SVGPathSeg.PATHSEG_UNKNOWN
  }, i.prototype.initialCommandIsMoveTo = function() {
   if (!this.hasMoreData()) return !0;
   var t = this.peekSegmentType();
   return t == window.SVGPathSeg.PATHSEG_MOVETO_ABS || t == window.SVGPathSeg.PATHSEG_MOVETO_REL
  }, i.prototype._parseNumber = function() {
   var t = 0,
    e = 0,
    i = 1,
    n = 0,
    r = 1,
    o = 1,
    s = this._currentIndex;
   if (this._skipOptionalSpaces(), this._currentIndex < this._endIndex && "+" == this._string.charAt(this._currentIndex) ? this._currentIndex++ : this._currentIndex < this._endIndex && "-" == this._string.charAt(this._currentIndex) && (this._currentIndex++, r = -1), !(this._currentIndex == this._endIndex || (this._string.charAt(this._currentIndex) < "0" || "9" < this._string.charAt(this._currentIndex)) && "." != this._string.charAt(this._currentIndex))) {
    for (var a = this._currentIndex; this._currentIndex < this._endIndex && "0" <= this._string.charAt(this._currentIndex) && this._string.charAt(this._currentIndex) <= "9";) this._currentIndex++;
    if (this._currentIndex != a)
     for (var h = this._currentIndex - 1, u = 1; a <= h;) e += u * (this._string.charAt(h--) - "0"), u *= 10;
    if (this._currentIndex < this._endIndex && "." == this._string.charAt(this._currentIndex)) {
     if (this._currentIndex++, this._currentIndex >= this._endIndex || this._string.charAt(this._currentIndex) < "0" || "9" < this._string.charAt(this._currentIndex)) return;
     for (; this._currentIndex < this._endIndex && "0" <= this._string.charAt(this._currentIndex) && this._string.charAt(this._currentIndex) <= "9";) i *= 10, n += (this._string.charAt(this._currentIndex) - "0") / i, this._currentIndex += 1
    }
    if (this._currentIndex != s && this._currentIndex + 1 < this._endIndex && ("e" == this._string.charAt(this._currentIndex) || "E" == this._string.charAt(this._currentIndex)) && "x" != this._string.charAt(this._currentIndex + 1) && "m" != this._string.charAt(this._currentIndex + 1)) {
     if (this._currentIndex++, "+" == this._string.charAt(this._currentIndex) ? this._currentIndex++ : "-" == this._string.charAt(this._currentIndex) && (this._currentIndex++, o = -1), this._currentIndex >= this._endIndex || this._string.charAt(this._currentIndex) < "0" || "9" < this._string.charAt(this._currentIndex)) return;
     for (; this._currentIndex < this._endIndex && "0" <= this._string.charAt(this._currentIndex) && this._string.charAt(this._currentIndex) <= "9";) t *= 10, t += this._string.charAt(this._currentIndex) - "0", this._currentIndex++
    }
    var c = e + n;
    if (c *= r, t && (c *= Math.pow(10, o * t)), s != this._currentIndex) return this._skipOptionalSpacesOrDelimiter(), c
   }
  }, i.prototype._parseArcFlag = function() {
   if (!(this._currentIndex >= this._endIndex)) {
    var t = !1,
     e = this._string.charAt(this._currentIndex++);
    if ("0" == e) t = !1;
    else {
     if ("1" != e) return;
     t = !0
    }
    return this._skipOptionalSpacesOrDelimiter(), t
   }
  }, i.prototype.parseSegment = function() {
   var t = this._string[this._currentIndex],
    e = this._pathSegTypeFromChar(t);
   if (e == window.SVGPathSeg.PATHSEG_UNKNOWN) {
    if (this._previousCommand == window.SVGPathSeg.PATHSEG_UNKNOWN) return null;
    if ((e = this._nextCommandHelper(t, this._previousCommand)) == window.SVGPathSeg.PATHSEG_UNKNOWN) return null
   } else this._currentIndex++;
   switch (this._previousCommand = e) {
    case window.SVGPathSeg.PATHSEG_MOVETO_REL:
     return new window.SVGPathSegMovetoRel(n, this._parseNumber(), this._parseNumber());
    case window.SVGPathSeg.PATHSEG_MOVETO_ABS:
     return new window.SVGPathSegMovetoAbs(n, this._parseNumber(), this._parseNumber());
    case window.SVGPathSeg.PATHSEG_LINETO_REL:
     return new window.SVGPathSegLinetoRel(n, this._parseNumber(), this._parseNumber());
    case window.SVGPathSeg.PATHSEG_LINETO_ABS:
     return new window.SVGPathSegLinetoAbs(n, this._parseNumber(), this._parseNumber());
    case window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_REL:
     return new window.SVGPathSegLinetoHorizontalRel(n, this._parseNumber());
    case window.SVGPathSeg.PATHSEG_LINETO_HORIZONTAL_ABS:
     return new window.SVGPathSegLinetoHorizontalAbs(n, this._parseNumber());
    case window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_REL:
     return new window.SVGPathSegLinetoVerticalRel(n, this._parseNumber());
    case window.SVGPathSeg.PATHSEG_LINETO_VERTICAL_ABS:
     return new window.SVGPathSegLinetoVerticalAbs(n, this._parseNumber());
    case window.SVGPathSeg.PATHSEG_CLOSEPATH:
     return this._skipOptionalSpaces(), new window.SVGPathSegClosePath(n);
    case window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_REL:
     var i = {
      x1: this._parseNumber(),
      y1: this._parseNumber(),
      x2: this._parseNumber(),
      y2: this._parseNumber(),
      x: this._parseNumber(),
      y: this._parseNumber()
     };
     return new window.SVGPathSegCurvetoCubicRel(n, i.x, i.y, i.x1, i.y1, i.x2, i.y2);
    case window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_ABS:
     return i = {
      x1: this._parseNumber(),
      y1: this._parseNumber(),
      x2: this._parseNumber(),
      y2: this._parseNumber(),
      x: this._parseNumber(),
      y: this._parseNumber()
     }, new window.SVGPathSegCurvetoCubicAbs(n, i.x, i.y, i.x1, i.y1, i.x2, i.y2);
    case window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_REL:
     return i = {
      x2: this._parseNumber(),
      y2: this._parseNumber(),
      x: this._parseNumber(),
      y: this._parseNumber()
     }, new window.SVGPathSegCurvetoCubicSmoothRel(n, i.x, i.y, i.x2, i.y2);
    case window.SVGPathSeg.PATHSEG_CURVETO_CUBIC_SMOOTH_ABS:
     return i = {
      x2: this._parseNumber(),
      y2: this._parseNumber(),
      x: this._parseNumber(),
      y: this._parseNumber()
     }, new window.SVGPathSegCurvetoCubicSmoothAbs(n, i.x, i.y, i.x2, i.y2);
    case window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_REL:
     return i = {
      x1: this._parseNumber(),
      y1: this._parseNumber(),
      x: this._parseNumber(),
      y: this._parseNumber()
     }, new window.SVGPathSegCurvetoQuadraticRel(n, i.x, i.y, i.x1, i.y1);
    case window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_ABS:
     return i = {
      x1: this._parseNumber(),
      y1: this._parseNumber(),
      x: this._parseNumber(),
      y: this._parseNumber()
     }, new window.SVGPathSegCurvetoQuadraticAbs(n, i.x, i.y, i.x1, i.y1);
    case window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_REL:
     return new window.SVGPathSegCurvetoQuadraticSmoothRel(n, this._parseNumber(), this._parseNumber());
    case window.SVGPathSeg.PATHSEG_CURVETO_QUADRATIC_SMOOTH_ABS:
     return new window.SVGPathSegCurvetoQuadraticSmoothAbs(n, this._parseNumber(), this._parseNumber());
    case window.SVGPathSeg.PATHSEG_ARC_REL:
     return i = {
      x1: this._parseNumber(),
      y1: this._parseNumber(),
      arcAngle: this._parseNumber(),
      arcLarge: this._parseArcFlag(),
      arcSweep: this._parseArcFlag(),
      x: this._parseNumber(),
      y: this._parseNumber()
     }, new window.SVGPathSegArcRel(n, i.x, i.y, i.x1, i.y1, i.arcAngle, i.arcLarge, i.arcSweep);
    case window.SVGPathSeg.PATHSEG_ARC_ABS:
     return i = {
      x1: this._parseNumber(),
      y1: this._parseNumber(),
      arcAngle: this._parseNumber(),
      arcLarge: this._parseArcFlag(),
      arcSweep: this._parseArcFlag(),
      x: this._parseNumber(),
      y: this._parseNumber()
     }, new window.SVGPathSegArcAbs(n, i.x, i.y, i.x1, i.y1, i.arcAngle, i.arcLarge, i.arcSweep);
    default:
     throw "Unknown path seg type."
   }
  };
  var r = new e,
   o = new i(t);
  if (!o.initialCommandIsMoveTo()) return [];
  for (; o.hasMoreData();) {
   var s = o.parseSegment();
   if (!s) return [];
   r.appendSegment(s)
  }
  return r.pathSegList
 })
}();
var __awaiter = this && this.__awaiter || function(t, s, a, h) {
  return new(a = a || Promise)(function(i, e) {
   function n(t) {
    try {
     o(h.next(t))
    } catch (t) {
     e(t)
    }
   }

   function r(t) {
    try {
     o(h.throw(t))
    } catch (t) {
     e(t)
    }
   }

   function o(t) {
    var e;
    t.done ? i(t.value) : ((e = t.value) instanceof a ? e : new a(function(t) {
     t(e)
    })).then(n, r)
   }
   o((h = h.apply(t, s || [])).next())
  })
 },
 __generator = this && this.__generator || function(i, n) {
  var r, o, s, t, a = {
   label: 0,
   sent: function() {
    if (1 & s[0]) throw s[1];
    return s[1]
   },
   trys: [],
   ops: []
  };
  return t = {
   next: e(0),
   throw: e(1),
   return: e(2)
  }, "function" == typeof Symbol && (t[Symbol.iterator] = function() {
   return this
  }), t;

  function e(e) {
   return function(t) {
    return function(e) {
     if (r) throw new TypeError("Generator is already executing.");
     for (; a;) try {
      if (r = 1, o && (s = 2 & e[0] ? o.return : e[0] ? o.throw || ((s = o.return) && s.call(o), 0) : o.next) && !(s = s.call(o, e[1])).done) return s;
      switch (o = 0, s && (e = [2 & e[0], s.value]), e[0]) {
       case 0:
       case 1:
        s = e;
        break;
       case 4:
        return a.label++, {
         value: e[1],
         done: !1
        };
       case 5:
        a.label++, o = e[1], e = [0];
        continue;
       case 7:
        e = a.ops.pop(), a.trys.pop();
        continue;
       default:
        if (!(s = 0 < (s = a.trys).length && s[s.length - 1]) && (6 === e[0] || 2 === e[0])) {
         a = 0;
         continue
        }
        if (3 === e[0] && (!s || e[1] > s[0] && e[1] < s[3])) {
         a.label = e[1];
         break
        }
        if (6 === e[0] && a.label < s[1]) {
         a.label = s[1], s = e;
         break
        }
        if (s && a.label < s[2]) {
         a.label = s[2], a.ops.push(e);
         break
        }
        s[2] && a.ops.pop(), a.trys.pop();
        continue
      }
      e = n.call(i, a)
     } catch (t) {
      e = [6, t], o = 0
     } finally {
      r = s = 0
     }
     if (5 & e[0]) throw e[1];
     return {
      value: e[0] ? e[1] : void 0,
      done: !0
     }
    }([e, t])
   }
  }
 },
 PdbTopologyViewerPlugin = function() {
  function t() {
   this.defaultColours = {
    domainSelection: "rgb(255,0,0)",
    mouseOver: "rgb(105,105,105)",
    borderColor: "rgb(0,0,0)",
    qualityGreen: "rgb(0,182.85714285714286,0)",
    qualityRed: "rgb(291.42857142857144,0,0)",
    qualityYellow: "rgb(364.2857142857143,364.2857142857143,75.71428571428572)",
    qualityRiboVision: "rgb(364.2857142857143,364.2857142857143,75.71428571428572)",
    qualityOrange: "rgb(291.42857142857144,121.42857142857143,0)"
   }, this.displayStyle = "border:1px solid #696969;", this.errorStyle = "border:1px solid #696969; height:54%; padding-top:46%; text-align:center; font-weight:bold;", this.menuStyle = "position:relative;height:38px;line-height:38px;background-color:#696969;padding: 0 10px;font-size:16px; color: #efefef;", this.svgWidth = 100, this.svgHeight = 100, this.subscribeEvents = !0, this.createNewEvent = function(t) {
    var n = {};
    return t.forEach(function(t, e) {
     var i;
     "function" == typeof MouseEvent ? i = new MouseEvent(t, {
      view: window,
      bubbles: !0,
      cancelable: !0
     }) : "function" == typeof document.createEvent && (i = document.createEvent("MouseEvents")).initEvent(t, !0, !0), n[t] = i
    }), n
   }, this.getAnnotationFromMappings = function() {
    var r = this;
    if (void 0 !== this.apiData[1]) {
     function t(t) {
      if (void 0 !== o[s[t]] && 0 !== Object.entries(o[s[t]]).length) {
       var e = [],
        i = o[s[t]];
       for (var n in i) i[n].mappings.forEach(function(t) {
        t.entity_id == r.entityId && t.chain_id == r.chainId && t.entropy_id == r.entropyId && e.push({
         start: t.start.residue_number,
         end: t.end.residue_number,
         color: void 0
        })
       });
       0 < e.length && a.domainTypes.push({
        label: s[t],
        data: e
       })
      }
     }
     for (var o = this.apiData[1][this.entryId], s = ["UniProt", "CATH", "Pfam", "SCOP", "RiboVision", "RiboVision1"], a = this, e = 0; e < 3; e++) t(e)
    }
   }, this.createDomainDropdown = function() {
    if (this.domainTypes = [{
      label: "Annotation",
      data: null
     }], this.getAnnotationFromMappings(), this.getAnnotationFromOutliers(),this.getAnnotationFromRiboVision(),this.getAnnotationFromRiboVision1(),this.getAnnotationFromRiboVision2(),this.getAnnotationFromRiboVision3(), this.selectedDomain = this.domainTypes[0], 1 < this.domainTypes.length) {
     var i = "";
	 this.domainTypes.forEach(function(t, e) {
      i = i + '<option value="' + e + '">' + t.label + "</option>"
     });
     var t = this.targetEle.querySelector(".menuSelectbox");
     t.innerHTML = i, t.addEventListener("change", this.displayDomain.bind(this)), this.targetEle.querySelector(".resetIcon").addEventListener("click", this.resetDisplay.bind(this)), this.targetEle.querySelector(".saveSVG").addEventListener("click", this.saveSVG.bind(this))
    } else this.targetEle.querySelector(".menuOptions").style.display = "none"
   }
  }
  return t.prototype.render = function(t, e) {
   var i = this;
   //Added this here
   this.entropyId = t.entropyId
   e && void 0 !== e.displayStyle && null != e.displayStyle && (this.displayStyle += e.displayStyle), e && void 0 !== e.errorStyle && null != e.errorStyle && (this.errorStyle += e.errorStyle), e && void 0 !== e.menuStyle && null != e.menuStyle && (this.menuStyle += e.menuStyle), this.targetEle = t, this.targetEle && (this.targetEle.innerHTML = ""), t && e && e.entryId && e.entityId ? (0 == e.subscribeEvents && (this.subscribeEvents = !1), this.entityId = e.entityId, this.entryId = e.entryId.toLowerCase(), void 0 === e.chainId || null == e.chainId ? this.getObservedResidues(this.entryId).then(function(t) {
    void 0 !== t && void 0 !== t[i.entryId] && void 0 !== t[i.entryId][i.entityId] ? (i.chainId = t[i.entryId][i.entityId][0].chain_id, i.initPainting()) : i.displayError()
   }) : (this.chainId = e.chainId, this.initPainting())) : this.displayError("param")
  }, t.prototype.initPainting = function() {
   var e = this;
   this.getApiData(this.entryId, this.chainId, this.entropyId).then(function(t) {
    if (t) {
     if (void 0 === t[0] || void 0 === t[2] || void 0 === t[4] || void 0 === t[5]) return void e.displayError();
     e.apiData = t, e.pdbevents = e.createNewEvent(["PDB.topologyViewer.click", "PDB.topologyViewer.mouseover", "PDB.topologyViewer.mouseout"]), e.getPDBSequenceArray(e.apiData[0][e.entryId]), e.drawTopologyStructures(), e.createDomainDropdown(), e.subscribeEvents && e.subscribeWcEvents()
    }
   })
  }, t.prototype.displayError = function(t) {
   var e = "Error: Data not available!";
   "param" == t && (e = "Error: Invalid Parameters!"), this.targetEle && (this.targetEle.innerHTML = '<div style="' + this.errorStyle + '">' + e + "</div>")
  }, t.prototype.getObservedResidues = function(i) {
   return __awaiter(this, void 0, void 0, function() {
    var e;
    return __generator(this, function(t) {
     switch (t.label) {
      case 0:
       return t.trys.push([0, 3, , 4]), [4, fetch("https://www.ebi.ac.uk/pdbe/api/pdb/entry/observed_residues_ratio/" + i)];
      case 1:
       return [4, t.sent().json()];
      case 2:
       return [2, t.sent()];
      case 3:
       return e = t.sent(), console.log("Couldn't load UniProt variants", e), [3, 4];
      case 4:
       return [2]
     }
    })
   })
  }, t.prototype.getApiData = function(i, n, entr) {
   return __awaiter(this, void 0, void 0, function() {
    var e;
    return __generator(this, function(t) {
     return e = ["https://www.ebi.ac.uk/pdbe/api/pdb/entry/entities/" + i, "https://www.ebi.ac.uk/pdbe/api/mappings/" + i, "https://www.ebi.ac.uk/pdbe/api/topology/entry/" + i, "http://127.0.0.1:8000/static/alignments/Data_dummy.txt", "https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/" + i + "/chain/" + n, "http://127.0.0.1:8000/alignments/" + entr,], [2, Promise.all(e.map(function(t) {
      return fetch(t)
     })).then(function(t) {
      return Promise.all(t.map(function(t) {
       return 200 == t.status ? t.json() : void 0
      }))
     })]
    })
   })
  }, t.prototype.getPDBSequenceArray = function(t) {
   for (var e = t.length, i = 0; i < e; i++) t[i].entity_id == this.entityId && (this.sequenceArr = t[i].sequence.split(""))
  }, t.prototype.chunkArray = function(t, e) {
   for (var i = [], n = 0, r = t.length; n < r;) i.push(t.slice(n, n += e));
   return i
  }, t.prototype.getDomainRange = function() {
   var e = this,
    i = [],
    t = this.apiData[2][this.entryId][this.entityId][this.chainId];
   for (var n in t) t[n] && t[n].forEach(function(t) {
    void 0 !== t.path && 0 < t.path.length && (i = i.concat(e.chunkArray(t.path, 2)))
   });
   this.xScale = d3.scaleLinear().domain([d3.min(i, function(t) {
    return t[0]
   }), d3.max(i, function(t) {
    return t[0]
   })]).range([1, this.svgWidth - 1]), this.yScale = d3.scaleLinear().domain([d3.min(i, function(t) {
    return t[1]
   }), d3.max(i, function(t) {
    return t[1]
   })]).range([1, this.svgHeight - 1]), this.zoom = d3.zoom().on("zoom", function() {
    return e.zoomDraw()
   })
  }, t.prototype.drawStrandSubpaths = function(t, e, i) {
   for (var n = this, r = e - t + 1, o = (this.scaledPointsArr[7] - this.scaledPointsArr[1]) / r, s = [], a = 0; a < r; a++) {
    var h = {
     type: "strands",
     elementIndex: i
    };
    0 === a ? (h.residue_number = t, h.pathData = [this.scaledPointsArr[4], this.scaledPointsArr[1], this.scaledPointsArr[4], this.scaledPointsArr[1] + o, this.scaledPointsArr[8], this.scaledPointsArr[1] + o, this.scaledPointsArr[8], this.scaledPointsArr[13]]) : (h.residue_number = t + a, h.pathData = [s[a - 1].pathData[2], s[a - 1].pathData[3], s[a - 1].pathData[2], s[a - 1].pathData[3] + o, s[a - 1].pathData[4], s[a - 1].pathData[5] + o, s[a - 1].pathData[4], s[a - 1].pathData[5]]), s.push(h)
   }
   this.svgEle.selectAll(".subpath-strands" + i).remove(), this.svgEle.selectAll(".subpath-strands" + i).data(s).enter().append("path").attr("class", function(t, e) {
    //console.log("strandsSubPath subpath-strands" + i + " topo_res_" + t.residue_number);
	return "strandsSubPath subpath-strands" + i + " topo_res_" + t.residue_number
   }).attr("d", function(t, e) {
    return "M " + t.pathData.join(" ") + " Z"
   }).attr("stroke", "#111").attr("stroke-width", "0").attr("fill", "white").attr("fill-opacity", "0").on("mouseover", function(t) {
    n.mouseoverAction(this, t)
   }).on("mousemove", function(t) {
    n.mouseoverAction(this, t)
   }).on("mouseout", function(t) {
    n.mouseoutAction(this, t)
   }).on("click", function(t) {
    n.clickAction(t)
   })
  }, t.prototype.drawStrandMaskShape = function(i) {
   var n = this.scaledPointsArr,
    t = [7, 8, 10, 12],
    e = [0, 1, 2, 3, 4, 5, 9, 11, 13];
   n[0] > n[6] && (t = [0, 1, 2, 3, 4, 5, 9, 11, 13], e = [7, 8, 10, 12]);
   for (var r = t.length, o = 0; o < r; o++) n[t[o]] = n[t[o]] + .3;
   for (var s = e.length, a = 0; a < s; a++) n[e[a]] = n[e[a]] - .3;
   n[14] = n[8], n[15] = n[13], n[16] = n[8], n[17] = n[7], n[18] = n[4], n[19] = n[7], n[20] = n[4], n[21] = n[1], this.svgEle.selectAll(".maskpath-strands" + i).remove(), this.svgEle.selectAll(".maskpath-strands" + i).data([n]).enter().append("path").attr("class", function(t, e) {
    return "strandMaskPath maskpath-strands" + i
   }).attr("d", function(t, e) {
    return "M" + n.join(" ") + "Z"
   }).attr("stroke", "#111").attr("stroke-width", .3).attr("fill", "white").attr("stroke-opacity", 0)
  }, t.prototype.renderTooltip = function(t, e) {
   var i = d3.select(".pdbTopologyTooltip");
   if (null == i._groups[0][0] && (i = d3.select("body").append("div").attr("class", "pdbTopologyTooltip").attr("style", "display: none;width: auto;position: absolute;background: #fff;padding: 5px;border: 1px solid #666;border-radius: 5px;box-shadow: 5px 6px 5px 0 rgba(0,0,0,.17);font-size: .9em;color: #555;z-index: 998;")), "show" === e) {
    var n = d3.event.pageX,
     r = d3.event.pageY,
     o = "Residue " + t.residue_number + " (" + this.sequenceArr[t.residue_number - 1] + ")";
    void 0 !== t.tooltipMsg && (o = void 0 !== t.tooltipPosition && "postfix" === t.tooltipPosition ? o + " " + t.tooltipMsg : t.tooltipMsg + " " + o), i.html(o).style("display", "block").style("top", r + 15 + "px").style("left", n + 10 + "px")
   } else i.style("display", "none")
  }, t.prototype.dispatchEvent = function(t, e, i) {
   var n = this.targetEle;
   void 0 !== i && (n = i), void 0 !== e && (this.pdbevents[t].eventData = e), n.dispatchEvent(this.pdbevents[t])
  }, t.prototype.clickAction = function(t) {
   this.dispatchEvent("PDB.topologyViewer.click", {
    residueNumber: t.residue_number,
    type: t.type,
    entropyId: this.entropyId,
    entryId: this.entryId,
    entityId: this.entityId,
    chainId: this.chainId
   })
  }, t.prototype.mouseoverAction = function(t, e) {
   var i = d3.select(t);
   this.renderTooltip(e, "show"), "strands" !== e.type && "helices" !== e.type || i.attr("fill", this.defaultColours.mouseOver).attr("fill-opacity", "0.3"), "coils" === e.type && i.attr("stroke", this.defaultColours.mouseOver).attr("stroke-width", 1.0), this.dispatchEvent("PDB.topologyViewer.mouseover", {
    residueNumber: e.residue_number,
    type: e.type,
    entropyId: this.entropyId,
    entryId: this.entryId,
    entityId: this.entityId,
    chainId: this.chainId
   })
  }, t.prototype.mouseoutAction = function(t, e) {
   var i = "white",
    n = 0,
    r = .3,
    o = d3.select(t);
   this.renderTooltip("", "hide"), o.classed("coloured") ? (i = o.attr("data-color"), r = n = 1) : "coils" === e.type && (i = this.defaultColours.borderColor), "strands" !== e.type && "helices" !== e.type || o.attr("fill", i).attr("fill-opacity", n), "coils" === e.type && o.attr("stroke", i).attr("stroke-width", r), this.dispatchEvent("PDB.topologyViewer.mouseout", {
    entropyId: this.entropyId,
    entryId: this.entryId,
    entityId: this.entityId,
    chainId: this.chainId
   })
  }, t.prototype.drawHelicesSubpaths = function(t, e, i, n) {
   var r = this,
    o = -5;
   this.scaledPointsArr[3] > this.scaledPointsArr[9] && (o = 5);
   var s = e - t + 1;
   o = 0;
   var a = (this.scaledPointsArr[9] - o - this.scaledPointsArr[3]) / s,
    h = 0,
    u = this.svgEle.select(".helices" + i).node().getBBox().height + a / 2,
    c = u / s;
   h = (a = (u - c / 2) / s) - c / 10, this.scaledPointsArr[3] > this.scaledPointsArr[9] && (h = -(u + c));
   for (var d = [], l = {}, p = 0; p < s; p++) l = {
    type: "helices"
   }, 0 === p ? (this.scaledPointsArr[3] > this.scaledPointsArr[9] ? l.residue_number = e : l.residue_number = t, l.pathData = [this.scaledPointsArr[0], this.scaledPointsArr[3] + h, this.scaledPointsArr[4], this.scaledPointsArr[3] + h, this.scaledPointsArr[4], this.scaledPointsArr[3] + h + a, this.scaledPointsArr[0], this.scaledPointsArr[3] + h + a]) : (this.scaledPointsArr[3] > this.scaledPointsArr[9] ? l.residue_number = e - p : l.residue_number = t + p, l.pathData = [d[p - 1].pathData[6], d[p - 1].pathData[7], d[p - 1].pathData[4], d[p - 1].pathData[5], d[p - 1].pathData[4], d[p - 1].pathData[5] + a, d[p - 1].pathData[6], d[p - 1].pathData[5] + a]), d.push(l);
   this.svgEle.selectAll(".subpath-helices" + i).remove(), this.svgEle.selectAll(".subpath-helices" + i).data(d).enter().append("path").attr("class", function(t) {
    return "helicesSubPath subpath-helices" + i + " topo_res_" + t.residue_number
   }).attr("d", function(t) {
    return "M" + t.pathData.join(" ") + " Z"
   }).attr("stroke", "#111").attr("stroke-width", "0").attr("fill", "white").attr("fill-opacity", "0").on("mouseover", function(t) {
    r.mouseoverAction(this, t)
   }).on("mousemove", function(t) {
    r.mouseoverAction(this, t)
   }).on("mouseout", function(t) {
    r.mouseoutAction(this, t)
   }).on("click", function(t) {
    r.clickAction(t)
   })
  }, t.prototype.drawHelicesMaskShape = function(e) {
   var t = [
    [this.scaledPointsArr[0] - .3, this.scaledPointsArr[1], this.scaledPointsArr[2], this.scaledPointsArr[3] - .3, this.scaledPointsArr[4] + .3, this.scaledPointsArr[5], this.scaledPointsArr[4] + .3, this.scaledPointsArr[3], this.scaledPointsArr[0] - .3, this.scaledPointsArr[3]],
    [this.scaledPointsArr[6] + .3, this.scaledPointsArr[7], this.scaledPointsArr[8], this.scaledPointsArr[9] + .3, this.scaledPointsArr[10] - .3, this.scaledPointsArr[11], this.scaledPointsArr[10] - .3, this.scaledPointsArr[9], this.scaledPointsArr[6] + .3, this.scaledPointsArr[9]]
   ];
   this.scaledPointsArr[3] > this.scaledPointsArr[9] && (t = [
    [this.scaledPointsArr[0] - .3, this.scaledPointsArr[1], this.scaledPointsArr[2], this.scaledPointsArr[3] + 2, this.scaledPointsArr[4] + .3, this.scaledPointsArr[5], this.scaledPointsArr[4] + .3, this.scaledPointsArr[3], this.scaledPointsArr[0] - .3, this.scaledPointsArr[3]],
    [this.scaledPointsArr[6] + .3, this.scaledPointsArr[7], this.scaledPointsArr[8], this.scaledPointsArr[9] - .3, this.scaledPointsArr[10] - .3, this.scaledPointsArr[11], this.scaledPointsArr[10] - .3, this.scaledPointsArr[9], this.scaledPointsArr[6] + .3, this.scaledPointsArr[9]]
   ]), this.svgEle.selectAll(".maskpath-helices" + e).remove(), this.svgEle.selectAll(".maskpath-helices" + e).data(t).enter().append("path").attr("class", function(t) {
    return "helicesMaskPath maskpath-helices" + e
   }).attr("d", function(t) {
    return "M" + t[0] + " " + t[1] + " Q" + t[2] + " " + t[3] + " " + t[4] + " " + t[5] + " L" + t[6] + " " + t[7] + " " + t[8] + " " + t[9] + " Z"
   }).attr("stroke", "#111").attr("stroke-width", .3).attr("fill", "white").attr("stroke-opacity", 0)
  }, t.prototype.drawCoilsSubpaths = function(t, e, i) {
   var n = this,
    r = this.svgEle.select(".coils" + i),
    o = e - t + 1,
    s = r.node().getTotalLength() / o,
    a = [],
    h = void 0,
    u = void 0,
    c = {};
   if (1 == o) c = {
    residue_number: t,
    type: "coils",
    pathData: n.scaledPointsArr,
    elementIndex: i
   }, a.push(c);
   else
    for (var d = 0; d < o; d++) {
     var l = s * (d + 1),
      p = r.node().getPointAtLength(l),
      S = r.node().getPathSegAtLength(l);
     c = {
      residue_number: t + d,
      type: "coils",
      elementIndex: i
     }, 1 === S ? c.pathData = n.scaledPointsArr.slice(0, 2) : (void 0 === u ? c.pathData = n.scaledPointsArr.slice(0, 2 * S) : (c.pathData = n.scaledPointsArr.slice(2 * u, 2 * S), c.pathData.unshift(h.x, h.y)), h = p, u = S), c.pathData = c.pathData.concat([p.x, p.y]), a.push(c)
    } - 1 !== t && -1 !== e && (this.svgEle.selectAll(".subpath-coils" + i).remove(), this.svgEle.selectAll(".subpath-coils" + i).data(a).enter().append("path").attr("class", function(t) {
     return "coilsSubPath subpath-coils" + i + " topo_res_" + t.residue_number
    }).attr("d", function(t) {
     return "M " + t.pathData.join(" ")
    }).attr("stroke", this.defaultColours.borderColor).attr("stroke-width", .3).attr("fill", "none").attr("stroke-opacity", "1").on("mouseover", function(t) {
     n.mouseoverAction(this, t)
    }).on("mousemove", function(t) {
     n.mouseoverAction(this, t)
    }).on("mouseout", function(t) {
     n.mouseoutAction(this, t)
    }).on("click", function(t) {
     n.clickAction(t)
    }), this.svgEle.selectAll(".coils" + i).attr("stroke-opacity", 0));
   var g = this.apiData[2][this.entryId][this.entityId][this.chainId].terms,
    w = this.apiData[2][this.entryId][this.entityId][this.chainId].coils.length;
   if (0 === i) this.svgEle.selectAll(".terminal_N").remove(), this.svgEle.selectAll(".terminal_N").data([g[0]]).enter().append("text").attr("class", "terminals terminal_N").attr("text-anchor", "middle").text("N").attr("x", a[0].pathData[0]).attr("y", a[0].pathData[1]).attr("stroke", "#0000ff").attr("stroke-width", "0.3").attr("font-size", "3px").attr("style", "-webkit-tap-highlight-color: rgba(0, 0, 0, 0); text-anchor: middle; font-style: normal; font-variant: normal; font-weight: normal; font-stretch: normal; line-height: normal; font-family: Arial;");
   else if (i === w - 1) {
    var _ = a[o - 1].pathData.length,
     y = -2;
    a[o - 1].pathData[_ - 1] > a[o - 1].pathData[_ - 3] && (y = 2), this.svgEle.selectAll(".terminal_C").remove(), this.svgEle.selectAll(".terminal_C").data([g[1]]).enter().append("text").attr("class", "terminals terminal_N").attr("text-anchor", "middle").text("C").attr("x", a[o - 1].pathData[_ - 2]).attr("y", a[o - 1].pathData[_ - 1] + y).attr("stroke", "#ff0000").attr("stroke-width", "0.3").attr("font-size", "3px").attr("style", "-webkit-tap-highlight-color: rgba(0, 0, 0, 0); text-anchor: middle; font-style: normal; font-variant: normal; font-weight: normal; font-stretch: normal; line-height: normal; font-family: Arial;")
   }
  }, t.prototype.drawTopologyStructures = function() {
   var u = this;
   this.targetEle.innerHTML = '<div style="' + this.displayStyle + '">\n            <div class="svgSection" style="position:relative;width:100%;"></div>\n            <div style="' + this.menuStyle + '">\n               <img src="https://www.ebi.ac.uk/pdbe/entry/static/images/logos/PDBe/logo_T_64.png" style="height:15px; width: 15px; border:0;position: absolute;margin-top: 11px;" />\n                <a style="color: #efefef;border-bottom:none; cursor:pointer;margin-left: 16px;" target="_blank" href="https://pdbe.org/' + this.entryId + '">' + this.entryId + '</a> | <span class="menuDesc">Entity ' + this.entityId + " | Chain " + this.chainId.toUpperCase() + '</span>\n                <div class="menuOptions" style="float:right;margin-right: 20px;">\n                    <select class="menuSelectbox" style="margin-right: 20px;"><option value="">Select</option></select>\n                <img class="saveSVG" src="http://apollo2.chemistry.gatech.edu/RiboVision3/pdb-topology-viewer-master_2/build/Save.png" style="height:15px; width: 15px; border:0;position: relative;margin-right: 15px;cursor:pointer;" title="saveSVG" />\n    <img class="resetIcon" src="https://www.ebi.ac.uk/pdbe/pdb-component-library/images/refresh.png" style="height:15px; width: 15px; border:0;position: absolute;margin-top: 11px;cursor:pointer;" title="Reset view" />\n             </div>\n            </div>\n        </div>';
   
   
   var t = this.targetEle.offsetWidth,
    e = this.targetEle.offsetHeight;
   0 == t && (t = this.targetEle.parentElement.offsetWidth), 0 == e && (e = this.targetEle.parentElement.offsetHeight), t <= 330 && (this.targetEle.querySelector(".menuDesc").innerText = this.entityId + " | " + this.chainId.toUpperCase());
   var i = this.targetEle.querySelector(".svgSection"),
    n = e - 40,
    r = t;
   i.style.height = n + "px";
   var o = n - 20,
    s = r - 5;

   function a(h) {
    var t = c[h];
    if (!t) return {
     value: void 0
    };
    t.forEach(function(a, t) {
     if (void 0 !== a.path && 0 < a.path.length && "terms" !== h) {
      var e = 0;
      if ("helices" === h) {
       var i = a.path[0] + (a.path[2] - a.path[0]) / 2;
       e = 1.3 * a.minoraxis * 2, a.path[1] > a.path[3] && (e = 1.3 * a.minoraxis * -2);
       var n = [a.path[0], a.path[1], i, a.path[1] - e, a.path[2], a.path[1], a.path[2], a.path[3], i, a.path[3] + e, a.path[0], a.path[3]];
       a.path = n
      }
      a.secStrType = h, a.pathIndex = t;
      var r = u.svgEle.selectAll("path." + h + t).data([a]).enter().append("path").attr("class", function() {
       return -1 === a.start && -1 === a.stop && "terms" !== h ? "dashedEle topologyEle " + h + " " + h + t + " topoEleRange_" + a.start + "-" + a.stop : "topologyEle " + h + " " + h + t + " topoEleRange_" + a.start + "-" + a.stop
      }).attr("d", function(t) {
       for (var e = "M", i = a.path.length, n = !0, r = 0; r < i; r++) {
        if ("helices" !== h || 2 !== r && 8 !== r || (e += " Q"), ("helices" === h && 6 === r || "coils" === h && a.path.length < 12 && 8 === r) && (e += " L"), n) {
         var o = u.xScale(a.path[r]);
         e += " " + o, u.scaledPointsArr.push(o)
        } else {
         var s = u.yScale(a.path[r]);
         e += " " + s, u.scaledPointsArr.push(s)
        }
        n = !n
       }
       return "strands" !== h && "helices" !== h || (e += " Z"), e
      }).attr("fill", "none").attr("stroke-width", .3).attr("stroke", u.defaultColours.borderColor); - 1 === a.start && -1 === a.stop && r.attr("stroke-dasharray", "0.9"), "strands" === h && (u.drawStrandSubpaths(a.start, a.stop, t), u.drawStrandMaskShape(t), u.svgEle._groups[0][0].append(r.node())), "helices" === h && (u.drawHelicesSubpaths(a.start, a.stop, t, e), u.drawHelicesMaskShape(t), u.svgEle._groups[0][0].append(r.node())), "coils" === h && u.drawCoilsSubpaths(a.start, a.stop, t), u.scaledPointsArr = []
     }
    })
   }
   i.innerHTML = '<svg class="topoSvg" preserveAspectRatio="xMidYMid meet" viewBox="0 0 100 100" style="width:' + s + "px;height:" + o + 'px;margin:10px 0;"></svg>', this.svgEle = d3.select(this.targetEle).select(".topoSvg"), this.getDomainRange(), this.scaledPointsArr = [], this.svgEle.call(this.zoom).on("contextmenu", function(t, e) {
    d3.event.preventDefault()
   });
   var c = this.apiData[2][this.entryId][this.entityId][this.chainId];
   for (var h in c) {
    var d = a(h);
    if ("object" == typeof d) return d.value
   }
   this.svgEle._groups[0][0].append(this.svgEle.selectAll(".validationResidue").node())
  }, t.prototype.zoomDraw = function() {
   var e = this,
    a = d3.event.transform.rescaleX(this.xScale),
    h = d3.event.transform.rescaleY(this.yScale),
    u = this;
   u.scaledPointsArr = [];
   var t = this.svgEle.selectAll(".topologyEle"),
    c = 0,
    d = 0,
    l = 0;
   t.each(function(t) {
    d3.select(d3.select(this).node()).attr("d", function(t) {
     c = t.pathIndex, d = t.start, l = t.stop;
     for (var e = "M", i = t.path.length, n = !0, r = 0; r < i; r++) {
      if ("helices" !== t.secStrType || 2 !== r && 8 !== r || (e += " Q"), ("helices" === t.secStrType && 6 === r || "coils" === t.secStrType && t.path.length < 12 && 8 === r) && (e += " L"), n) {
       var o = a(t.path[r]);
       e += " " + o, u.scaledPointsArr.push(o)
      } else {
       var s = h(t.path[r]);
       e += " " + s, u.scaledPointsArr.push(s)
      }
      n = !n
     }
     return "strands" !== t.secStrType && "helices" !== t.secStrType || (e += " Z"), e
    }), "helices" === t.secStrType ? (u.drawHelicesSubpaths(d, l, c, 0), u.drawHelicesMaskShape(c), u.svgEle._groups[0][0].append(d3.select(this).node())) : "strands" === t.secStrType ? (u.drawStrandSubpaths(d, l, c), u.drawStrandMaskShape(c), u.svgEle._groups[0][0].append(d3.select(this).node())) : "coils" === t.secStrType && u.drawCoilsSubpaths(d, l, c), u.scaledPointsArr = []
   });
   var s = 0;
   this.svgEle.selectAll(".validationResidue").attr("transform", function(t) {
    var e = u.svgEle.select(".topo_res_" + t.residue_number),
     i = e.node().getBBox(),
     n = e.data(),
     r = {
      x: 0,
      y: 0
     };
    if ("strands" === n[0].type || "helices" === n[0].type) r = {
     x: i.x + i.width / 2,
     y: i.y + i.height / 2
    };
    else {
     var o = e.node().getPointAtLength(e.node().getTotalLength() / 2);
     r = {
      x: o.x,
      y: o.y
     }
    }
    return s = i.height / 2, "translate(" + r.x + "," + r.y + ")"
   }).attr("d", d3.symbol().type(function(t, e) {
    return d3.symbols[0]
   }).size(s)), this.svgEle.selectAll(".residueSelection").attr("d", function(t) {
    var e = d3.select(this).data();
    return u.svgEle.select(".topo_res_" + e[0].residueNumber).attr("d")
   }), this.svgEle._groups[0][0].querySelectorAll(".coilsSubPath").forEach(function(t) {
    return e.svgEle._groups[0][0].append(t)
   }), this.svgEle._groups[0][0].querySelectorAll(".dashedEle").forEach(function(t) {
    return e.svgEle._groups[0][0].append(t)
   }), this.displayDomain("zoom"), this.svgEle._groups[0][0].querySelectorAll(".validationResidue").forEach(function(t) {
    return e.svgEle._groups[0][0].append(t)
   }), this.svgEle._groups[0][0].querySelectorAll(".residueSelection").forEach(function(t) {
    return e.svgEle._groups[0][0].append(t)
   })
  }, t.prototype.clearHighlight = function() {
   this.svgEle.selectAll(".residueHighlight").remove()
  }, t.prototype.highlight = function(t, e, r, o) {
   function i(e) {
    var t = d.svgEle.select(".topo_res_" + e);
    if (t && t._groups && null == t._groups[0][0]) return {
     value: void 0
    };
    var i = t.node(),
     n = t.data();
    r && (a = "string" == typeof r ? h = r : (h = d3.rgb(r.r, r.g, r.b), d3.rgb(r.r, r.g, r.b))), "strands" !== n[0].type && "helices" !== n[0].type ? (a = "none", u = 2, c = .5) : h = "none", d.svgEle.append("path").data([{
     residueNumber: e
    }]).attr("class", function(t) {
     return "click" == o ? "residueSelection seletectedResidue_" + e : "residueHighlight highlightResidue_" + e
    }).attr("d", t.attr("d")).attr("fill", a).attr("fill-opacity", .5).attr("stroke", h).attr("stroke-opacity", c).attr("stroke-width", u).on("mouseover", function(t) {
     s.mouseoverAction(i, n[0])
    }).on("mousemove", function(t) {
     s.mouseoverAction(i, n[0])
    }).on("mouseout", function(t) {
     s.mouseoutAction(i, n[0])
    }).on("click", function(t) {
     s.clickAction(n[0])
    })
   }
   for (var s = this, a = "#000000", h = "#000000", u = .3, c = 0, d = this, n = t; n <= e; n++) {
    var l = i(n);
    if ("object" == typeof l) return l.value
   }
  }, t.prototype.drawValidationShape = function(t, e, i) {
   var n = this,
    r = n.svgEle.select(".topo_res_" + t);
   if (null != r._groups[0][0]) {
    var o = r.node().getBBox(),
     s = r.data(),
     a = {
      x: 0,
      y: 0
     };
    if ("strands" === s[0].type || "helices" === s[0].type) a = {
     x: o.x + o.width / 2,
     y: o.y + o.height / 2
    };
    else {
     var h = r.node().getPointAtLength(r.node().getTotalLength() / 2);
     a = {
      x: h.x,
      y: h.y
     }
    }
    var u = {
     residue_number: t,
     tooltipMsg: "Validation issue: RSRZ <br>",
     tooltipPosition: "prefix"
    };
    this.svgEle.append("path").attr("class", "validationResidue rsrz_" + t).data([u]).attr("fill", i).attr("stroke", "#000").attr("stroke-width", .1).attr("transform", function(t) {
     return "translate(" + a.x + "," + a.y + ")"
    }).attr("d", d3.symbol().type(function(t, e) {
     return d3.symbols[0]
    }).size(1.8)).style("display", "none").on("mouseover", function(t) {
     n.mouseoverAction(this, t)
    }).on("mousemove", function(t) {
     n.mouseoverAction(this, t)
    }).on("mouseout", function(t) {
     n.mouseoutAction(this, t)
    }).on("click", function(t) {
     n.clickAction(t)
    })
   }
  }, t.prototype.getChainStartAndEnd = function() {
   if (void 0 !== this.apiData[4]) {
    for (var t = this.apiData[4][this.entryId].molecules[0].chains, i = {
      start: 0,
      end: 0
     }, e = t.length, n = 0; n < e; n++)
     if (t[n].chain_id == this.chainId) {
      t[n].observed.forEach(function(t, e) {
       0 == e ? (i.start = t.start.residue_number, i.end = t.end.residue_number) : (t.start.residue_number < i.start && (i.start = t.start.residue_number), t.end.residue_number > i.end && (i.end = t.end.residue_number))
      });
      break
     } return i
   }
  }, t.prototype.getAnnotationFromOutliers = function() {
   var e = this,
    o = this,
    t = this.getChainStartAndEnd(),

    s = [{
     start: t.start,
     end: t.end,
  
     tooltipMsg: "No validation issue reported for "
    }],
    a = [],
    h = [0];
	console.log(t);
   if (void 0 !== this.apiData[3]) {
    var i = this.apiData[3][this.entryId];
    void 0 !== i && void 0 !== i.molecules && 0 < i.molecules.length && (i.molecules.forEach(function(t) {
     t.entity_id == e.entityId && t.chains.forEach(function(t) {
      t.chain_id == e.chainId && t.models.forEach(function(t) {
       t.residues.forEach(function(t) {
        var e = o.defaultColours.qualityYellow,
         i = "issue2";
        if (1 !== t.outlier_types.length || "RSRZ" !== t.outlier_types[0]) {
         1 === t.outlier_types.length ? e = o.defaultColours.qualityOrange : (e = 2 === t.outlier_types.length ? o.defaultColours.qualityOrange : o.defaultColours.qualityRed, i = "issues"), h.push(t.residue_number);
         var n = "Validation " + i + ": " + t.outlier_types.join(", ") + "<br>"; - 1 < a.indexOf(t.residue_number) && (n = "Validation issues: " + t.outlier_types.join(", ") + ", RSRZ<br>"), s.push({
          start: parseInt(t.residue_number),
          end: parseInt(t.residue_number),
          
          tooltipMsg: n,
          tooltipPosition: "prefix"
         })
        } else {
         o.drawValidationShape(t.residue_number, "circle", e), a.push(t.residue_number);
         var r = h.indexOf(t.residue_number); - 1 < r ? s[r].tooltipMsg = s[r].tooltipMsg.replace("<br>", ", RSRZ<br>") : (s.push({
          start: parseInt(t.residue_number),
          end: parseInt(t.residue_number),
        
          tooltipMsg: "Validation issue: RSRZ <br>",
          tooltipPosition: "prefix"
         }), h.push(t.residue_number))
        }
       })
      })
     })
    }), 0 < s.length && this.domainTypes.push({
     label: "Quality",
     data: s
    }))
   }
  },
  
  t.prototype.getAnnotationFromRiboVision = function() {
    var e = this,
     o = this,
     t = this.getChainStartAndEnd(),
 
     s = [{
      start: t.start,
      end: t.end,
      color: o.defaultColours.qualityGreen,
      tooltipMsg: "No validation issue reported for "
     }],
     a = [],
     h = [0];
     console.log(t);
     var RiboData=this.apiData[5];
     //console.log(RiboData);
    if (void 0 !== this.apiData[5]) {
     var i = this.apiData[3]["1b23"];
     //console.log(this.apiData[3], this.entryId);
     var RiboData=this.apiData[5];
     //console.log(this.apiData[5], this.apiData[5].length,  i.molecules);


     RiboData_Y=getCol(this.apiData[5], 1);
   

     Y_min=Math.min.apply(Math, RiboData_Y);
     Y_max=Math.max.apply(Math, RiboData_Y);
     Y_range=Math.max.apply(Math, RiboData_Y)-Math.min.apply(Math, RiboData_Y);

     console.log(Y_range);

     RiboData_Y_norm = RiboData_Y.map(x => (x-Y_min)/Y_range);
    

 
     void 0 !== i && void 0 !== i.molecules && 0 < i.molecules.length && (i.molecules.forEach(function(t) {
      t.entity_id == e.entityId && t.chains.forEach(function(t) {
       t.chain_id == e.chainId && t.models.forEach(function(t) {
         t.residues.forEach(function(t) {
         console.log(t.residue_number); 
         if (RiboData.length > t.residue_number) {console.log(t.residue_number, RiboData[parseInt(t.residue_number)][1]); 
            
            var col_tol=palette.tolDivergingColor(RiboData_Y_norm[parseInt(t.residue_number)]);

            var col_tol_rgb= hexToRgb1(col_tol);
           
         
         //o.defaultColours.qualityRiboVision= "rgb(364.2857142857143,364.2857142857143,75.71428571428572)"
         o.defaultColours.qualityRiboVision= "rgb("+String(col_tol_rgb)+")"

        // o.defaultColours.qualityRiboVision= "rgb("+String(RiboData[parseInt(t.residue_number)][1]*100)+","+String(RiboData[parseInt(t.residue_number)][1]*100)+",175.71428571428572)"
         var e = o.defaultColours.qualityRiboVision,

    

          i = "issue";
          console.log(e)
         }
         if ((1 !== t.outlier_types.length || "RSRZ" !== t.outlier_types[0] ) && RiboData.length > t.residue_number) {
          1 === t.outlier_types.length ? e = "rgb("+String(col_tol_rgb)+")" : (e = 2 === t.outlier_types.length ? o.defaultColours.qualityOrange : o.defaultColours.qualityRed, i = "issues"), h.push(t.residue_number);
          var n = "Validation " + i + ": " + t.outlier_types.join(", ") + "<br>"; - 1 < a.indexOf(t.residue_number) && (n = "Validation issues: " + t.outlier_types.join(", ") + ", RSRZ<br>"), s.push({
           start: parseInt(t.residue_number),
           end: parseInt(t.residue_number),
           color: e,
           tooltipMsg: n,
           tooltipPosition: "prefix"
          })
         } else {
          e = o.defaultColours.qualityRed, o.drawValidationShape(t.residue_number, "circle", e), a.push(t.residue_number);
          var r = h.indexOf(t.residue_number); - 1 < r ? s[r].tooltipMsg = s[r].tooltipMsg.replace("<br>", ", RSRZ<br>") : (s.push({
           start: parseInt(t.residue_number),
           end: parseInt(t.residue_number),
           color: o.defaultColours.qualityGreen,
           tooltipMsg: "Validation issue: RSRZ <br>",
           tooltipPosition: "prefix"
          }), h.push(t.residue_number))
         }
        })
       })
      })
     }), 0 < s.length && this.domainTypes.push({
      label: "RiboVision",
      data: s
     }))
    }
   },

   t.prototype.getAnnotationFromRiboVision2 = function() {
    var e = this,
     o = this,
     t = this.getChainStartAndEnd(),
 
     s = [{
      start: t.start,
      end: t.end,
   
      tooltipMsg: "No validation issue reported for "
     }],
     a = [],
     h = [0];
     console.log(t);
    if (void 0 !== this.apiData[3]) {
     var i = this.apiData[3]["1b23"];
     void 0 !== i && void 0 !== i.molecules && 0 < i.molecules.length && (i.molecules.forEach(function(t) {
      t.entity_id == e.entityId && t.chains.forEach(function(t) {
       t.chain_id == e.chainId && t.models.forEach(function(t) {
        t.residues.forEach(function(t) {
         var e = o.defaultColours.qualityYellow,
          i = "issue2";
         if (1 !== t.outlier_types.length || "RSRZ" !== t.outlier_types[0]) {
          1 === t.outlier_types.length ? e = o.defaultColours.qualityOrange : (e = 2 === t.outlier_types.length ? o.defaultColours.qualityOrange : o.defaultColours.qualityRed, i = "issues"), h.push(t.residue_number);
          var n = "Validation " + i + ": " + t.outlier_types.join(", ") + "<br>"; - 1 < a.indexOf(t.residue_number) && (n = "Validation issues: " + t.outlier_types.join(", ") + ", RSRZ<br>"), s.push({
           start: parseInt(t.residue_number),
           end: parseInt(t.residue_number),
           
           tooltipMsg: n,
           tooltipPosition: "prefix"
          })
         } else {
          o.drawValidationShape(t.residue_number, "circle", e), a.push(t.residue_number);
          var r = h.indexOf(t.residue_number); - 1 < r ? s[r].tooltipMsg = s[r].tooltipMsg.replace("<br>", ", RSRZ<br>") : (s.push({
           start: parseInt(t.residue_number),
           end: parseInt(t.residue_number),
         
           tooltipMsg: "Validation issue: RSRZ <br>",
           tooltipPosition: "prefix"
          }), h.push(t.residue_number))
         }
        })
       })
      })
     }), 0 < s.length && this.domainTypes.push({
      label: "RiboVision2",
      data: s
     }))
    }
   },
   t.prototype.getAnnotationFromRiboVision1 = function() {
    var e = this,
     o = this,
     t = this.getChainStartAndEnd(),
 
     s = [{
      start: t.start,
      end: t.end,
      
      tooltipMsg: "No validation issue reported for "
     }],
     a = [],
     h = [0];
     console.log(t);
     var RiboData=this.apiData[5];
    
    if (void 0 !== this.apiData[5]) {
     var i = this.apiData[3]["1b23"];
     var RiboData=this.apiData[5];

     RiboData_Y=getCol(this.apiData[5], 1);
   

     Y_min=Math.min.apply(Math, RiboData_Y);
     Y_max=Math.max.apply(Math, RiboData_Y);
     Y_range=Math.max.apply(Math, RiboData_Y)-Math.min.apply(Math, RiboData_Y);

     console.log(Y_range);

     RiboData_Y_norm = RiboData_Y.map(x => (x-Y_min)/Y_range);
    




   


  
 
     void 0 !== i && void 0 !== i.molecules && 0 < i.molecules.length && (i.molecules.forEach(function(t) {
      t.entity_id == "2" && t.chains.forEach(function(t) {
       t.chain_id == "P" && t.models.forEach(function(t) {
         t.residues.forEach(function(t) {
        
         if (RiboData.length > t.residue_number) {
             //console.log(t.residue_number, RiboData[parseInt(t.residue_number)][1]); 

             RiboData[parseInt(t.residue_number)][1]


             var col_tol=palette.tolDivergingColor(RiboData_Y_norm[parseInt(t.residue_number)]);

             var col_tol_rgb= hexToRgb1(col_tol);
             //console.log(col_tol, col_tol_rgb);
            
            
         
         //o.defaultColours.qualityRiboVision= "rgb(364.2857142857143,364.2857142857143,75.71428571428572)"
        // o.defaultColours.qualityRiboVision= "rgb("+String(RiboData[parseInt(t.residue_number)][1]*100)+","+String(RiboData[parseInt(t.residue_number)][1]*100)+",175.71428571428572)"

         o.defaultColours.qualityRiboVision= "rgb("+String(col_tol_rgb)+")"
         var e = o.defaultColours.qualityRiboVision,

    

          i = "issue";
          //console.log(e, t.outlier_types.length)
         }
         if ((1 === t.outlier_types.length  ) && RiboData.length > t.residue_number) {
          //e = "rgb("+String(RiboData[parseInt(t.residue_number)][1]*10)+","+String(RiboData[parseInt(t.residue_number)][1]*100)+",175.71428571428572)", o.drawValidationShape(t.residue_number, "circle", o.defaultColours.qualityRiboVision); 
          e = "rgb("+String(col_tol_rgb)+")", o.drawValidationShape(t.residue_number, "circle", o.defaultColours.qualityRiboVision); 
          var l = i.indexOf(t.residue_number);
          //console.log(l);

         
          if (l ===-1) {s.push({
            start: parseInt(t.residue_number),
            end: parseInt(t.residue_number),
            color: e,
            tooltipMsg: Number.parseFloat(RiboData[parseInt(t.residue_number)][1]).toPrecision(3),
            tooltipPosition: "prefix"
          }),  h.push(t.residue_number);
          o.drawValidationShape(t.residue_number, "circle", e);
      
          a.push(t.residue_number)
          }
          
         } 
        })
       })
      })
     }), 0 < s.length && this.domainTypes.push({
      label: "RiboVision1",
      data: s
     }))
    }
   }, 

   t.prototype.getAnnotationFromRiboVision3 = function() {
    var e = this,
     o = this,
     t = this.getChainStartAndEnd(),
 
     s = [{
      start: t.start,
      end: t.end,
      
      tooltipMsg: "No validation issue reported for "
     }],
     a = [],
     h = [0];
     console.log(t);
     var RiboData=this.apiData[5];
    
    if (void 0 !== this.apiData[5]) {
     var i = this.apiData[3]["1b23"];
     var RiboData=this.apiData[5];

     RiboData_Y=getCol(this.apiData[5], 1);
   

     Y_min=Math.min.apply(Math, RiboData_Y);
     Y_max=Math.max.apply(Math, RiboData_Y);
     Y_range=Math.max.apply(Math, RiboData_Y)-Math.min.apply(Math, RiboData_Y);

     console.log(Y_range);

     RiboData_Y_norm = RiboData_Y.map(x => (x-Y_min)/Y_range);
    




   


  
 
     void 0 !== i && void 0 !== i.molecules && 0 < i.molecules.length && (i.molecules.forEach(function(t) {
      t.entity_id == "2" && t.chains.forEach(function(t) {
       t.chain_id == "P" && t.models.forEach(function(t) {
         t.residues.forEach(function(t) {
        
         if (RiboData.length > t.residue_number) {
             //console.log(t.residue_number, RiboData[parseInt(t.residue_number)][1]); 

             RiboData[parseInt(t.residue_number)][1]


             var col_tol=palette.tolDivergingColor(RiboData_Y_norm[parseInt(t.residue_number)]);

             var col_tol_rgb= hexToRgb1(col_tol);

             var col_tol_rgb_GB = hexToRgb_GreenBlind2(col_tol);
             console.log(col_tol, col_tol_rgb, col_tol_rgb_GB );

          
            
            
         
         //o.defaultColours.qualityRiboVision= "rgb(364.2857142857143,364.2857142857143,75.71428571428572)"
        // o.defaultColours.qualityRiboVision= "rgb("+String(RiboData[parseInt(t.residue_number)][1]*100)+","+String(RiboData[parseInt(t.residue_number)][1]*100)+",175.71428571428572)"

         o.defaultColours.qualityRiboVision= "rgb("+String(col_tol_rgb_GB)+")"
         var e = o.defaultColours.qualityRiboVision,

    

          i = "issue";
          //console.log(e, t.outlier_types.length)
         }
         if ((1 === t.outlier_types.length  ) && RiboData.length > t.residue_number) {
          //e = "rgb("+String(RiboData[parseInt(t.residue_number)][1]*10)+","+String(RiboData[parseInt(t.residue_number)][1]*100)+",175.71428571428572)", o.drawValidationShape(t.residue_number, "circle", o.defaultColours.qualityRiboVision); 
          e = "rgb("+String(col_tol_rgb_GB)+")", o.drawValidationShape(t.residue_number, "circle", o.defaultColours.qualityRiboVision); 
          var l = i.indexOf(t.residue_number);
          //console.log(l);

         
          if (l ===-1) {s.push({
            start: parseInt(t.residue_number),
            end: parseInt(t.residue_number),
            color: e,
            tooltipMsg: "Validation issue: RSRZ1 <br>",
            tooltipPosition: "prefix"
          }),  h.push(t.residue_number);
          o.drawValidationShape(t.residue_number, "circle", e);
      
          a.push(t.residue_number)
          }
          
         } 
        })
       })
      })
     }), 0 < s.length && this.domainTypes.push({
      label: "RiboVision3",
      data: s
     }))
    }
   },  
  t.prototype.resetTheme = function() {
   var o = this;
   this.svgEle.selectAll(".coloured").each(function(t) {
    var e = d3.select(this),
     i = e.node();
    e.data()[0].tooltipMsg = void 0, e.data()[0].tooltipPosition = void 0;
    var n = d3.select(i).classed("coloured", !1),
     r = n.attr("class").split(" "); - 1 < r.indexOf("strandsSubPath") || -1 < r.indexOf("helicesSubPath") ? n.attr("fill", "white").attr("fill-opacity", 0) : n.attr("stroke", o.defaultColours.borderColor).attr("stroke-width", .3)
   }), this.svgEle.selectAll(".validationResidue").style("display", "none")
  }, t.prototype.changeResidueColor = function(t, e, i, n) {
   void 0 === e && (e = this.defaultColours.domainSelection);
   var r = this.svgEle.select(".topo_res_" + t);
   null != r._groups[0][0] && (r.data()[0].tooltipMsg = i, r.data()[0].tooltipPosition = n, r.attr("stroke", function(t) {
    return "coils" === t.type ? e : "#111"
   }).attr("stroke-width", function(t) {
    return "coils" === t.type ? 0.4 : 0
   }).attr("fill", function(t) {
    return "coils" === t.type ? "none" : e
   }).attr("fill-opacity", function(t) {
    return "coils" === t.type ? 0 : 1
   }).classed("coloured", !0).attr("data-color", e))
  }, t.prototype.updateTheme = function(t) {
   var i = this;
   t.forEach(function(t) {
    for (var e = t.start; e <= t.end; e++) i.changeResidueColor(e, t.color, t.tooltipMsg, t.tooltipPosition)
   })
  }, 
  t.prototype.updateTheme_RV2 = function(t) {
    var i = this;
    t.forEach(function(t) {
     for (var e = t.start; e <= 1; e++) i.changeResidueColor(e, t.color, t.tooltipMsg, t.tooltipPosition)
    })
   }, 
  
  t.prototype.displayDomain = function(t) {
   var e = this.targetEle.querySelector(".menuSelectbox"),
    i = parseInt(e.value),
    n = this.domainTypes[i];
    console.log(n.data);
    null !== n.data ? (this.resetTheme(), this.updateTheme(n.data), "Quality" === n.label && this.svgEle.selectAll(".validationResidue").style("display", "block")) : "zoom" !== t && this.resetTheme();
    if ("RiboVision2" === n.label) {null !== n.data ? (this.resetTheme(), this.updateTheme_RV2(n.data), "RiboVision2" === n.label && this.svgEle.selectAll(".validationResidue").style("display", "block")) : "zoom" !== t && this.resetTheme()} 
   

  },  
  
  t.prototype.resetDisplay = function() {
   this.targetEle.querySelector(".menuSelectbox").value = 0, this.displayDomain()
  }, t.prototype.saveSVG = function() {
        function getNode(n, v) {
		  n = document.createElementNS("http://www.w3.org/2000/svg", n);
		   for (var p in v) 
		   n.setAttributeNS(null, p.replace(/[0-9]/g,'o').replace(/\$/g,'d').replace(/\[/g,'b').replace(/[A-Z]/g, function(m, p, o, s) { return "-" + m.toLowerCase(); }), v[p]);
		  return n
		}

		

		  var svgData1=this.targetEle.querySelector(".topoSvg")
	   

		


	    var svg = getNode("svg");
        document.body.appendChild(svg);
        
		svg.appendChild(svgData1);
		
		
		function saveSvg1(svgEl, name) {
         svgEl.setAttribute("xmlns", "http://www.w3.org/2000/svg");
         var svgData = svgEl.outerHTML;
         var preface = '<?xml version="1.0" standalone="no"?>\r\n';
          var svgBlob = new Blob([preface, svgData], {type:"image/svg+xml;charset=utf-8"});
         var svgUrl = URL.createObjectURL(svgBlob);
         var downloadLink = document.createElement("a");
          downloadLink.href = svgUrl;
          downloadLink.download = name;
          document.body.appendChild(downloadLink);
          downloadLink.click();
          document.body.removeChild(downloadLink);
        }
	    saveSvg1(svg, 'test.svg')	
  },t.prototype.handleSeqViewerEvents = function(t, e) {
   if (void 0 !== t.eventData) {
    if (t.eventData.entryId.toLowerCase() != this.entryId.toLowerCase() || t.eventData.entityId != this.entityId) return;
    if (t.eventData.elementData.pathData.chain_id && t.eventData.elementData.pathData.chain_id != this.chainId) return;
    var i = "residueSelection";
    "mouseover" == e && (i = "residueHighlight"), this.svgEle.selectAll("." + i).remove();
    var n = void 0,
     r = void 0;
    if (t.eventData.residueNumber ? (n = t.eventData.residueNumber, r = t.eventData.residueNumber) : t.eventData.elementData.pathData.start.residue_number && t.eventData.elementData.pathData.end.residue_number && (n = t.eventData.elementData.pathData.start.residue_number, r = t.eventData.elementData.pathData.end.residue_number), void 0 !== n && void 0 !== r) {
     var o;
     o = t.eventData.elementData.color && 1 == t.eventData.elementData.color.length ? t.eventData.elementData.color[0] : {
      r: t.eventData.elementData.color[0],
      g: t.eventData.elementData.color[1],
      b: t.eventData.elementData.color[2]
     }, this.highlight(n, r, o, e)
    }
   }
  }, t.prototype.handleProtvistaEvents = function(t, e) {
   if (void 0 !== t.detail) {
    var i = void 0,
     n = "residueSelection";
    if ("mouseover" == e && (n = "residueHighlight"), this.svgEle.selectAll("." + n).remove(), void 0 !== t.detail.feature) {
     if (void 0 !== t.detail.feature.accession) {
      var r = t.detail.feature.accession.split(" ");
      if ("Chain" == r[0] && r[1].toLowerCase() != this.chainId.toLowerCase()) return
     } - 1 < t.detail.trackIndex && t.detail.feature.locations && t.detail.feature.locations[0].fragments[t.detail.trackIndex].color && (i = t.detail.feature.locations[0].fragments[t.detail.trackIndex].color), void 0 === i && t.detail.feature.color && (i = t.detail.feature.color)
    }
    void 0 === i && t.detail.color && (i = t.detail.color), void 0 !== i && (i = /rgb/g.test(i) ? i.substring(4, i.length - 1).split(",") : [i]);
    var o = void 0;
    i && (o = 1 == i.length ? t.eventData.elementData.color[0] : {
     r: t.eventData.elementData.color[0],
     g: t.eventData.elementData.color[1],
     b: t.eventData.elementData.color[2]
    }), this.highlight(t.detail.start, t.detail.end, o, e)
   }
  }, t.prototype.handleMolstarEvents = function(t, e) {
   if (void 0 !== t.eventData && 0 < Object.keys(t.eventData).length) {
    var i = "residueSelection";
    if ("mouseover" == e && (i = "residueHighlight"), this.svgEle.selectAll("." + i).remove(), t.eventData.entry_id.toLowerCase() != this.entryId.toLowerCase() || t.eventData.entity_id != this.entityId) return;
    if (t.eventData.label_asym_id.toLowerCase() != this.chainId.toLowerCase()) return;
    this.highlight(t.eventData.seq_id, t.eventData.seq_id, void 0, e)
   }
  }, t.prototype.subscribeWcEvents = function() {
   var e = this;
   document.addEventListener("PDB.seqViewer.click", function(t) {
    e.handleSeqViewerEvents(t, "click")
   }), document.addEventListener("PDB.seqViewer.mouseover", function(t) {
    e.handleSeqViewerEvents(t, "mouseover")
   }), document.addEventListener("PDB.seqViewer.mouseout", function() {
    e.svgEle.selectAll(".residueHighlight").remove()
   }), document.addEventListener("PDB.litemol.click", function(t) {
    e.svgEle.selectAll(".residueSelection").remove(), t.eventData.entryId.toLowerCase() == e.entryId.toLowerCase() && t.eventData.entityId == e.entityId && t.eventData.chainId.toLowerCase() == e.chainId.toLowerCase() && e.highlight(t.eventData.residueNumber, t.eventData.residueNumber, void 0, "click")
   }), document.addEventListener("PDB.litemol.mouseover", function(t) {
    e.svgEle.selectAll(".residueHighlight").remove(), t.eventData.entryId.toLowerCase() == e.entryId.toLowerCase() && t.eventData.entityId == e.entityId && t.eventData.chainId.toLowerCase() == e.chainId.toLowerCase() && e.highlight(t.eventData.residueNumber, t.eventData.residueNumber, void 0, "mouseover")
   }), document.addEventListener("protvista-click", function(t) {
    e.handleProtvistaEvents(t, "click")
   }), document.addEventListener("protvista-mouseover", function(t) {
    e.handleProtvistaEvents(t, "mouseover")
   }), document.addEventListener("protvista-mouseout", function() {
    e.svgEle.selectAll(".residueHighlight").remove()
   }), document.addEventListener("PDB.molstar.click", function(t) {
    e.handleMolstarEvents(t, "click")
   }), document.addEventListener("PDB.molstar.mouseover", function(t) {
    e.handleMolstarEvents(t, "mouseover")
   }), document.addEventListener("PDB.molstar.mouseout", function() {
    e.svgEle.selectAll(".residueHighlight").remove()
   })
  }, t
 }();
window.PdbTopologyViewerPlugin = PdbTopologyViewerPlugin,
 function(i) {
  var n = {};

  function r(t) {
   if (n[t]) return n[t].exports;
   var e = n[t] = {
    i: t,
    l: !1,
    exports: {}
   };
   return i[t].call(e.exports, e, e.exports, r), e.l = !0, e.exports
  }
  r.m = i, r.c = n, r.d = function(t, e, i) {
   r.o(t, e) || Object.defineProperty(t, e, {
    enumerable: !0,
    get: i
   })
  }, r.r = function(t) {
   "undefined" != typeof Symbol && Symbol.toStringTag && Object.defineProperty(t, Symbol.toStringTag, {
    value: "Module"
   }), Object.defineProperty(t, "__esModule", {
    value: !0
   })
  }, r.t = function(e, t) {
   if (1 & t && (e = r(e)), 8 & t) return e;
   if (4 & t && "object" == typeof e && e && e.__esModule) return e;
   var i = Object.create(null);
   if (r.r(i), Object.defineProperty(i, "default", {
     enumerable: !0,
     value: e
    }), 2 & t && "string" != typeof e)
    for (var n in e) r.d(i, n, function(t) {
     return e[t]
    }.bind(null, n));
   return i
  }, r.n = function(t) {
   var e = t && t.__esModule ? function() {
    return t.default
   } : function() {
    return t
   };
   return r.d(e, "a", e), e
  }, r.o = function(t, e) {
   return Object.prototype.hasOwnProperty.call(t, e)
  }, r.p = "", r(r.s = 7)
 }([function(i, t) {
  function n(t, e) {
   return i.exports = n = Object.setPrototypeOf || function(t, e) {
    return t.__proto__ = e, t
   }, n(t, e)
  }
  i.exports = n
 }, function(e, t) {
  function i(t) {
   return e.exports = i = Object.setPrototypeOf ? Object.getPrototypeOf : function(t) {
    return t.__proto__ || Object.getPrototypeOf(t)
   }, i(t)
  }
  e.exports = i
 }, function(t, e) {
  function n(t, e) {
   for (var i = 0; i < e.length; i++) {
    var n = e[i];
    n.enumerable = n.enumerable || !1, n.configurable = !0, "value" in n && (n.writable = !0), Object.defineProperty(t, n.key, n)
   }
  }
  t.exports = function(t, e, i) {
   return e && n(t.prototype, e), i && n(t, i), t
  }
 }, function(t, e) {
  t.exports = function(t, e) {
   if (!(t instanceof e)) throw new TypeError("Cannot call a class as a function")
  }
 }, function(t, e, i) {
  var n = i(8),
   r = i(9);
  t.exports = function(t, e) {
   return !e || "object" !== n(e) && "function" != typeof e ? r(t) : e
  }
 }, function(t, e, i) {
  var n = i(0);
  t.exports = function(t, e) {
   if ("function" != typeof e && null !== e) throw new TypeError("Super expression must either be null or a function");
   t.prototype = Object.create(e && e.prototype, {
    constructor: {
     value: t,
     writable: !0,
     configurable: !0
    }
   }), e && n(t, e)
  }
 }, function(e, t, i) {
  var n = i(1),
   r = i(0),
   o = i(10),
   s = i(11);

  function a(t) {
   var i = "function" == typeof Map ? new Map : void 0;
   return e.exports = a = function(t) {
    if (null === t || !o(t)) return t;
    if ("function" != typeof t) throw new TypeError("Super expression must either be null or a function");
    if (void 0 !== i) {
     if (i.has(t)) return i.get(t);
     i.set(t, e)
    }

    function e() {
     return s(t, arguments, n(this).constructor)
    }
    return e.prototype = Object.create(t.prototype, {
     constructor: {
      value: e,
      enumerable: !1,
      writable: !0,
      configurable: !0
     }
    }), r(e, t)
   }, a(t)
  }
  e.exports = a
 }, function(t, e, i) {
  "use strict";
  i.r(e);
  var n, r = i(3),
   o = i.n(r),
   s = i(4),
   a = i.n(s),
   h = i(1),
   u = i.n(h),
   c = i(2),
   d = i.n(c),
   l = i(5),
   p = i.n(l),
   S = i(6),
   g = (n = i.n(S)()(HTMLElement), p()(w, n), d()(w, null, [{
    key: "observedAttributes",
    get: function() {
     return ["entry-id", "entity-id", "entropy-id", "chain-id", "display-style", "error-style", "menu-style", "subscribe-events"]
    }
   }]), d()(w, [{
    key: "validateParams",
    value: function() {
     return void 0 !== this.entryId && void 0 !== this.entityId && null != this.entityId
    }
   }, {
    key: "invokePlugin",
    value: function() {
     if (this.validateParams()) {
      void 0 === this.pluginInstance && (this.pluginInstance = new PdbTopologyViewerPlugin);
      var t = {
       entryId: this.entryId,
       entityId: this.entityId,
       entropyId: this.entropyId
      };
      void 0 !== this.entropyId && null !== this.entropyId && (t.entropyId = this.entropyId), void 0 !== this.chainId && null !== this.chainId && (t.chainId = this.chainId), void 0 !== this.displayStyle && null !== this.displayStyle && (t.displayStyle = this.displayStyle), void 0 !== this.errorStyle && null !== this.errorStyle && (t.errorStyle = this.errorStyle), void 0 !== this.menuStyle && null !== this.menuStyle && (t.menuStyle = this.menuStyle), void 0 !== this.subscribeEvents && null !== this.subscribeEvents && (t.subscribeEvents = this.subscribeEvents), this.pluginInstance.render(this, t)
     }
    }
   }, {
    key: "attributeChangedCallback",
    value: function() {
     this.entropyId = this.getAttribute("entropy-id"), this.entryId = this.getAttribute("entry-id"), this.entityId = this.getAttribute("entity-id"), this.chainId = this.getAttribute("chain-id"), this.displayStyle = this.getAttribute("display-style"), this.errorStyle = this.getAttribute("error-style"), this.menuStyle = this.getAttribute("menu-style"), this.subscribeEvents = this.getAttribute("subscribe-events"), this.invokePlugin()
    }
   }]), w);

  function w() {
   return o()(this, w), a()(this, u()(w).call(this))
  }
  e.default = g, customElements.define("pdb-topology-viewer", g)
 }, function(e, t) {
  function i(t) {
   return "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? e.exports = i = function(t) {
    return typeof t
   } : e.exports = i = function(t) {
    return t && "function" == typeof Symbol && t.constructor === Symbol && t !== Symbol.prototype ? "symbol" : typeof t
   }, i(t)
  }
  e.exports = i
 }, function(t, e) {
  t.exports = function(t) {
   if (void 0 === t) throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
   return t
  }
 }, function(t, e) {
  t.exports = function(t) {
   return -1 !== Function.toString.call(t).indexOf("[native code]")
  }
 }, function(n, t, e) {
  var o = e(0);

  function r(t, e, i) {
   return function() {
    if ("undefined" != typeof Reflect && Reflect.construct && !Reflect.construct.sham) {
     if ("function" == typeof Proxy) return 1;
     try {
      return Date.prototype.toString.call(Reflect.construct(Date, [], function() {})), 1
     } catch (t) {
      return
     }
    }
   }() ? n.exports = r = Reflect.construct : n.exports = r = function(t, e, i) {
    var n = [null];
    n.push.apply(n, e);
    var r = new(Function.bind.apply(t, n));
    return i && o(r, i.prototype), r
   }, r.apply(null, arguments)
  }
  n.exports = r
 }]);
