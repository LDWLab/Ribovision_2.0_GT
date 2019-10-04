$.widget( "ui.iosMenu", {
	options: {
		backText:      'Back',
		slideDuration: 200,
		slideEasing:   'linear'
	},

	_insertBackButtons: function() {
		this.element.find( 'li ul, li ol' ).prepend(
			$( '<li class="ui-menu-item ios-menu-back-link ui-corner-role" tabindex="-1" role="menuitem" aria-haspopup="true">' +
				'	<span class="ui-menu-icon ui-icon ui-icon-carat-1-w"></span>' +
							this.options.backText +
				 '</li>'
		) );
		return this;
	},

	_create: function( options ) {
		var iosMenu = this;

		iosMenu
			._insertBackButtons()
			.element
				.addClass( 'ios-style' )
				.menu({
					// When a submenu shows up, place it just to the right
					// of the current menu. Later, we'll slide it into view.
					position: {
						my: 'left top',
						at: 'right top',
						of: iosMenu.element
					}
				});

		//var menu = iosMenu.element.data( 'uiMenu' );
		var menu = $(iosMenu.element).data('ui-menu');
		// Override menu#select to account for nesting and back buttons:
		menu.select = function( event ) {
			//menu.active = menu.active || $( event.target ).closest( ".ui-menu-item" ); //new random line
			if ( menu.active && menu.active.hasClass("ios-menu-back-link") ) {
				// if you selected "back", go back:
				menu.focus( event, menu.active );
				menu.collapse( event );
				/*
				if ( menu.collapse( event ) ) {
					event.stopImmediatePropagation();
				}
				event.preventDefault();*/
			} else if ( menu.active && menu.active.find( '> ul' ).length ) {
				// if you selected something with children, show the children:
				menu.focus( event, menu.active );
				menu.expand( event );
				/*
				if ( menu.expand( event ) ) {
					event.stopImmediatePropagation();
				}
				event.preventDefault();*/
			} else {
				menu._trigger( 'select', event, { item: menu.active } );
			}
		};
		/*
		// Override menu#expand to add return true:
		menu.expand = function( event ) {
			var newItem = this.active &&
				this.active
					.children( ".ui-menu " )
					.children( ".ui-menu-item" )
					.first();

			if ( newItem && newItem.length ) {
				this._open( newItem.parent() );

				// Delay so Firefox will not hide activedescendant change in expanding submenu from AT
				this._delay(function() {
					this.focus( event, newItem );
				});
				return true;
			}
		};*/
	
		// Override menu#collapse to enable sliding behavior:
		menu.collapse = function( event ) {
			var newItem = this.active && this.active.parents( 'li:not(.ui-menubar-item) ').first(),
					self		= this,
					parent;
			if ( newItem && newItem.length ) {
			  newItem.addClass( 'ui-state-focus' ).removeClass( 'ui-state-active' );
				parent = this.active.parent();
				parent
					.attr( 'aria-hidden', 'true' )
					.attr( 'aria-expanded', 'false' )
					.animate({
						left: self.element.css( 'width' )
					}, iosMenu.options.slideDuration, iosMenu.options.slideEasing, function() {
						parent.hide();
						self.focus( event, newItem );
					})
				//return true;
			} else if ( event && event.which === $.ui.keyCode.ESCAPE ) {
				// #left gets called both for left-arrow and escape. If it's the
				// latter and we're at the top, fire a "close" event:
				self._trigger( 'close', event );
			}
		};

		
	
	
		// Override menu#_open to enable sliding behavior:
		var menuOpenWithoutSliding = menu._open;
		menu._open = function ( submenu ) {
			menuOpenWithoutSliding.call( this, submenu );

			submenu.animate({
				left: 0,
				height: menu.element[0].clientHeight - 4,
				width: menu.element[0].clientWidth 
			}, iosMenu.options.slideDuration, iosMenu.options.slideEasing);
		};

		// Override menu#_startOpening so that hovering doesn't
		// initiate the sliding:
		menu._startOpening = function( submenu ) {
			clearTimeout( this.timer );
			// Don't open if already open fixes a Firefox bug that caused a .5 pixel
			// shift in the submenu position when mousing over the carat icon
			if ( submenu.attr( "aria-hidden" ) !== "true" ) {
				return;
			}
		}
	},

	destroy: function() {
	  var menu = this.element && this.element.data( 'uiMenu' );
		menu && menu.destroy();
	}
});


