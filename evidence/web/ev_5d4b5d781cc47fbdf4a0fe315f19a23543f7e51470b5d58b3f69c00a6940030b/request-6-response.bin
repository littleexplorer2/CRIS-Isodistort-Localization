// isodistor.js
// Bob Hanson 2023.12.08
// BH 2024.01.28 removed html injection

appletCount = 0;

createIsoDistortApp = function(inputType, isoData, formData) {
	var app = (inputType == "isovizdistortion" ? "IsoDistortApp" : "IsoDiffractApp");
	var Info = {
			  code: null,
			  main: "org.byu.isodistort." + app,
			  args: [ isoData, formData, document ],
			  core: "_iso",//"NONE",//			
			  width: +isoData.appletwidth,
			  height: +isoData.appletwidth/1.6,
			  readyFunction: null,
			  serverURL: 'https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php',
			  j2sPath: 'swingjs/j2s',
			  console:'sysoutdiv',
			  allowjavascript: true
			}
	if (!$("#sysoutdiv")[0]) {
		$("div").last().append('<div style="width:600px;height:300px;">'
						+'<div id="sysoutdiv" contentEditable="true" style="border:1px solid green;width:100%;height:95%;overflow:auto"></div>'
						+'</div>');
	}
	var html = SwingJS.getAppletHtml("isoDistort" + (appletCount++), Info);	
	$($("p").get(1)).before(html);
}

fetchIsoDistort = function(inputType, callback) {
	// "isovizdistortion", createIsoDistortApplet
	// "isovizdiffraction", createIsoDistortApplet
	var formData = {};
	$("form input").each(function(i,e) {
		if (e.value == inputType) {
			e.checked = true;
		}
	});
	$("input").each(function(i,e) {
		if (e.name != 'appbtn' && (e.type != "radio" && e.type != "checkbox" || e.checked)) {
			formData[e.name] = e.value;			
		}
	});
	var f = $("form");
    var info = {
	      type: f.attr("method"),
	      url: f.attr("action"),
	      data: formData,
	      dataType: "text",
	      encode: true,
	    };
	$.ajax(info)
	      .done(function(data){callback(data, formData)})
	      .fail(function (data){alert("error submitting form: " + data)});
}

clickIsoDistort = function() {
	fetchIsoDistort("isovizdistortion", function(data, formData) { createIsoDistortApp("isovizdistortion", data, formData)});
}

