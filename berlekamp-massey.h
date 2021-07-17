#ifndef _BERLEKAMP_MASSEY_
#define _BERLEKAMP_MASSEY_
#include <vector>

namespace BM {
	const int mod = 1e9+7;

	long long add(long long a, long long b) {
		return ((a+b)%mod+mod)%mod;
	}

	long long powmod(long long a, long long b) {
		long long ret = 1;
		while(b) {
			if(b&1ll) ret*=a, ret%=mod;
			a*=a; a%=mod;
			b>>=1;
		}
		return ret;
	}

	long long inv(long long a) {
		return powmod(a, mod-2);
	}

	std::vector<long long> linearRecurrence(std::vector<long long> a) {
		int n = a.size();
		std::vector<long long> l = {0};
		std::vector<long long> mark = {0};
		std::vector<std::vector<long long> > co = {{1}};
		for(int i=1; i<n; i++) {
			long long eval = 0;
			for(int j=0; j<(int)l.size(); j++) {
				eval+=(l[j]*a[i-j-1]%mod+mod)%mod;
				eval%=mod;
			}
	 
			if(eval==a[i]) continue;

			//There is a discrepancy

			long long d = add(eval, -a[i]);
			int mn = n+2;
			int ind = -1;
			for(int j=0; j<(int)mark.size(); j++) {
				if(i-1-mark[j]+(int)co[j].size()<mn) {
					mn = i-1-mark[j]+(int)co[j].size();
					ind = j;
				}
			}


			mark.push_back(i);
			std::vector<long long> tmp = {1}; for(int x:l) tmp.push_back((mod-x)%mod); co.push_back(tmp);

			std::vector<long long> newl(i-1-mark[ind], 0); for(int x:co[ind]) newl.push_back(x);
			long long newval = 0; for(int j=0; j<(int)newl.size(); j++) newval = add(newval, (newl[j]*a[i-1-j]%mod+mod)%mod);
			for(long long &x:newl) x = x*inv(newval)%mod*d%mod;

			l.resize(std::max(l.size(), newl.size()));

			for(int j=0; j<(int)l.size(); j++) {
				if(j<(int)newl.size()) l[j] = add(l[j], -newl[j]);
			}		
		}
		return l;
	}

}

#endif