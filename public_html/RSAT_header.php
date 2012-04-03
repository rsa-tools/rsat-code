<html xmlns="http://www.w3.org/1999/xhtml">
	<head>
<?php
if ($result) {
	$res = " - results";
	$class = "results";
	echo '<link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">',"\n";
} else  {
	$res = "";
	$class = "form";
	echo '<link rel="stylesheet" type="text/css" href="main.css" media="screen,projection,print"/>',"\n";
}
?>
<!--  	<link rel="stylesheet" type="text/css" href="tabs.css" media="screen,projection,print"/>-->	
		<script src="RSAT_menu.js" type="text/javascript"></script>
		<script src="RSAT_tabs.js" type="text/javascript"></script>
		<title>RSAT - <?php echo $prog_name, $res; ?></title>
		<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	</head>

	<body class="<?php echo $class; ?>">
		<h3 align='center'><a href="<?php echo $properties['rsat_www']?>">RSA-tools</a> - <?php echo $prog_name, $res; ?></h3><br/>
