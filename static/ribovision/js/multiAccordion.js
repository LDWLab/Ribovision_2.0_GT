//Multi Accordion
(function($){

	$.fn.multiAccordion = function () {
		$(this).addClass("ui-accordion ui-accordion-icons ui-widget ui-helper-reset")
		.find("h3").not('.ui-accordion-header')
		.addClass("ui-accordion-header ui-helper-reset ui-state-default ui-corner-all ui-accordion-icons")
		
		.hover(function () {
			$(this).toggleClass("ui-state-hover");
		})
		.prepend('<span class="ui-accordion-header-icon ui-icon ui-icon-triangle-1-e"></span>')
		.click(function () {
			$(this)
			.toggleClass("ui-accordion-header-active ui-state-active ui-state-default ui-corner-bottom")
			.find("> .ui-icon").toggleClass("ui-icon-triangle-1-e ui-icon-triangle-1-s").end()
			.next().toggleClass("ui-accordion-content-active").slideToggle('fast'); //animation speed:fast
			return false;
		})	
		/*
		.dblclick(function(){
			alert("double clicked!");	
		})	*/
		.next()
		.addClass("ui-accordion-content ui-helper-reset ui-widget-content ui-corner-bottom")
		.css("display", "block")
		.hide()
		//.end().trigger("click");
	};

	$.fn.multiAccordionRefresh = function () {
		alert("biggie smalls");
	};
})(jQuery);