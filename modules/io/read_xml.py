def read_XML_script(xmlfile):
	
	"""
	Handler function to read the XML file and to create the GUI
	
	"""
	# global script_tree
	from lxml import etree
	import datetime
	import sys
	try:
		file_handle = open(xmlfile,'r')
	except FileNotFoundError:
		print(str(datetime.datetime.now())+" - "+xmlfile+': file not found.')
	else:
		try: 
			script_tree = etree.parse(xmlfile)
			#print("file opened... Je suis passé par là")
		except SyntaxError:
			sys.exit(": Invalid file...".format)
	
	return script_tree


