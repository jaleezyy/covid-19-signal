#!/usr/bin/env python

### Written by @fmaguire
### Modified for SIGNAL integration by @jaleezyy

import sys
import argparse
import subprocess
import urllib.request as web

"""
Example input version file would contain something like:

pangolin: v3.1.14
pangolearn: 2021-10-13
constellations: v0.0.18
scorpio: v0.3.13
pango-designation: v1.2.88
"""

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Update pangolin in conda env "
												 "to specific versions using "
												 "pip")
	parser.add_argument("--versions_file", required=True, type=str,
						help="File containing pangolin dependency "
									"versions e.g., \npangolin: 3.1.14"
									"\npangolearn: 2021-10-13"
									"\n...")
	args = parser.parse_args()

	# provides current pangolin install details
	# and load them in a dict for comparison
	print("## Existing pangolin install:")
	installed_versions = subprocess.run(["pangolin", "--all-versions"],
										check=True,
										stdout=subprocess.PIPE)
	installed_versions = installed_versions.stdout.decode('utf-8')
	print(installed_versions)

	installed_ver_dict = {}
	for dep_ver in map(str.strip, installed_versions.split('\n')):
		# skip empty line at end
		if len(dep_ver) == 0:
			continue

		try:
			dependency, version = dep_ver.split(': ')
		except ValueError:
			continue
			#print(dep_ver.split(': '))
			#dependency, version = dep_ver.split(': ')
		
		if dependency != "pangolearn" and not version.startswith("v"):
			version = "v" + version

		# tidy up pango-designation version for dep we actually change
		# instead of packages with pangolearn/usher models
		if dependency == 'pango-designation aliases':
			installed_ver_dict['pango-designation'] = version
		else:
			installed_ver_dict[dependency] = version

	# parse the dependency version file provided, validate real dependencies
	# tidy up version strings, and then use pip to update
	required = []
	valid_deps = ['pangolin', 'pangolearn', 'constellations',
				  'scorpio', 'pango-designation', 'pangolin-data']
	required_v3 = ['pangolin', 'scorpio', 'constellations', 'pangolearn', 'pango-designation'] 
	required_v4 = ['pangolin', 'scorpio', 'constellations', 'pangolin-data']
	print("## Changing installed versions as needed:")
	with open(args.versions_file) as fh:
		for line in fh:
			line = line.split(':')
			dependency = line[0].strip()
			requested_ver = line[1].strip()

			if dependency == "pangolin" and len(required) == 0:
				if requested_ver.startswith("3") or requested_ver.startswith("v3"):
					required = required_v3
				elif requested_ver.startswith("4") or requested_ver.startswith("v4") or str(requested_ver) == "None":
					required = required_v4
				else:
					required = valid_deps # no change to valid deps

			if dependency not in required:
				if dependency not in valid_deps:
					raise ValueError(f"{dependency} is not a valid pangolin "
									f"dependency. Must be in {required}")
				else:
					print(f"{dependency} not required! Skipping...")
					continue

			if dependency != "pangolearn" and not requested_ver.startswith("v"):
				requested_ver = "v" + requested_ver
			
			
			# Get latest version number to compare with installed when latest version is requested
			try:
				if (str(requested_ver) == "vNone") or (dependency == "pangolearn" and str(requested_ver) == "None"):
					commit_url = web.urlopen(f"https://github.com/cov-lineages/{dependency}/releases/latest").geturl()
					requested_ver = commit_url.split("/")[-1] # request version is latest
					try:
						installed_version = installed_ver_dict[dependency]
					except: # due to inconsistency, let's assume not installed
						installed_version = None
				else:
					try:
						installed_version = installed_ver_dict[dependency]
					except KeyError:
						installed_version = None
			except (web.HTTPError, web.URLError):
				print(f"Cannot determine latest version of {dependency}! Skipping update!")
				continue
			
			# if above block doesn't run, there is a specific version requested
			try:
				if requested_ver == installed_version:
					print(f"{dependency} not updated as requested {requested_ver} already installed")
				else:
					if installed_version is None:
						print(f"Installing {dependency} {requested_ver}")
					else:
						print(f"Changing {dependency} from {installed_version} to {requested_ver}")
					subprocess.run([sys.executable, '-m', 'pip', 'install',
								f"git+https://github.com/cov-lineages/{dependency}.git@{requested_ver}"],
								check=True,
								stdout=subprocess.DEVNULL,
								stderr=subprocess.DEVNULL)
			except (subprocess.CalledProcessError):
				print(f"Something went wrong updating {dependency}! Skipping update!")
				continue

	# provides pangolin install details after update to specific versions in supplied version file
	with open('final_pangolin_versions.txt', 'w+') as out:
		print("## Pangolin and dependencies now:", file=out)
		out_ver = subprocess.run(["pangolin", "--all-versions"], check=True, stdout=subprocess.PIPE)
		print(out_ver.stdout.decode("utf-8").replace("[32m****\nPangolin running in usher mode.\n****[0m", "").strip(), file=out)
