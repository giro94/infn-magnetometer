
vector<TString> getListOfFiles(TString folder){
	vector<TString> files;

	cout<<"Reading folder "<<folder<<"\n";
	DIR* dir = opendir(folder.Data());
	if (!dir) {
		cout<<"Cannot open folder "<<folder<<"!\n";
		throw std::runtime_error{"Cannot open folder"};
	}

	struct dirent* dirfile;
	while((dirfile = readdir(dir)) != NULL){
		if (dirfile->d_type != DT_REG) continue;
		TString fname = dirfile->d_name;
		if (!fname.EndsWith("csv")) continue;
		files.push_back(fname);
	}
	sort(files.begin(),files.end());

	return files;
}

TDatime getFileTime(TString filename){
	filename.ReplaceAll("FD_","");
	filename.ReplaceAll(".csv","");
	std::tm tm{};
	std::istringstream iss(filename.Data());
	iss >> std::get_time(&tm, "%m_%d_%Y %H_%M_%S");
	if (iss.fail()) {
		cout<<filename<<"\n";
		cout<<"Can't extract time from file name "<<filename<<"!\n";
	    throw std::runtime_error{"failed to parse time string"};
	}
	return TDatime(tm.tm_year+1900,tm.tm_mon,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
}


void plot_kicks(TString folder, TString output_file){


	//Begin reading of files
	vector<TString> files = getListOfFiles(folder);
	int Nfiles = files.size();

	
}