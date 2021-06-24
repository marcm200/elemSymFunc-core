struct DynSlowString {
	int32_t memory;
	char* text;

	DynSlowString(void);
	DynSlowString(const int32_t);
	virtual ~DynSlowString ();

	void add(const char*);
	void add(DynSlowString&);
	void add(const char);
	void setEmpty(void);
};

// DynSlowString
void DynSlowString::setEmpty(void) {
	if (text) text[0]=0;
}

DynSlowString::DynSlowString(const int32_t amem) {
	memory=amem;
	text=new char[memory];
	if (!text) {
		LOGMSG("\nError. Memory. DynSlowString.\n");
		exit(99);
	}
}

DynSlowString::DynSlowString() {
	memory=0;
	text=NULL;
}

DynSlowString::~DynSlowString() {
	if (text) delete[] text;
}

void DynSlowString::add(const char ac) {
	char tmp[16];
	sprintf(tmp,"%c",ac);

	add(tmp);
}

inline int32_t sum_int32t(
	const int32_t a,
	const int32_t b
) {
	int64_t sum=(int64_t)a + (int64_t)b;
	if (
		(sum < (-INT32MAX)) ||
		(sum >   INT32MAX)
	) {
		LOGMSG("\nError. Overflow by int32 additionﬂn");
		exit(99);
	}

	return (int32_t)sum; // lower half
}

void DynSlowString::add(const char* atext) {
	if (!atext) return;

	int32_t ttlen=strlen(atext);
	if (ttlen <= 0) return;

	int32_t currentlen;
	if (text) currentlen=strlen(text); else currentlen=0;

	if (
		(sum_int32t(currentlen,ttlen) > (memory-8) ) ||
		(!text)
	) {
		// new memory necessary
		int32_t mem0=sum_int32t(currentlen,ttlen);
		memory = sum_int32t(mem0,mem0);
		char* oldtext=text;
		text=new char[memory];
		if (!text) {
			LOGMSG("\nError. DynSlowString::add.\n");
			exit(99);
		}
		if (oldtext) {
			strcpy(text,oldtext);
			delete[] oldtext;
		} else {
			text[0]=0;
		}
	} // allocating

	// enough memory allocated
	sprintf(&text[currentlen],"%s",atext);
}

void DynSlowString::add(DynSlowString& as) {
	if (as.text) add(as.text);
}

