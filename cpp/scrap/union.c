
typedef union {
	double d;
	float f[2];
} pack;

int main() {
	pack p = {
		.d = 5
	}

}
