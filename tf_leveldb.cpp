#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<algorithm>
#include<vector>
#include<queue>
#include<cstring>
#include<string>
#include<map>
#include<cmath>
#include<stack>
#include<set>

using namespace std;
//ieeepaper
//最大编号：7102531 最小编号：15 总边数：10600587
//总结点数：22309488
const int maxN = 8e6;
const int maxn = 3e6;
const int maxm = 10700000;
int pre_to_new_NUM[maxN];
int tln[maxn]={0};//topological level number
int tot;
int tf[maxn*4];
int TF;
LeveldbGraph<Node, Edge> graph("mygraph.db");
vector<int> G[maxn];

vector<int>GG[maxn];
//int fa[maxn*4];
//***************************************tarjan算法求强连通分量********************************//
int pre[maxn]={0},lowlink[maxn],sccno[maxn],dfs_clock,scc_cnt;//sccno is the new number
void dfs(int u)  
{  
	static stack<int> s;
    pre[u] = lowlink[u] = ++dfs_clock;  
    s.push(u);  
    for(int i=0; i<GG[u].size(); i++)  
    {  
        int v=GG[u][i];  
        if(!pre[v])  
        {  
            dfs(v);  
            lowlink[u]=min(lowlink[u],lowlink[v]);  
        }  
        else if(!sccno[v])  
        {  
            lowlink[u]=min(lowlink[u],pre[v]);  
        }  
    }  
    if(lowlink[u]==pre[u])  
    {  
        scc_cnt++;  
        for(;;)  
        {  
            int x=s.top();  
            s.pop();  
            sccno[x]=scc_cnt;  
            if(x==u)  
                break;  
        }  
    }  
}  
void find_scc(int n)  
{  
    dfs_clock=scc_cnt=0;  
    memset(sccno,0,sizeof(sccno));   
    for(int i=0; i<n; i++)  
        if(!pre[i]) dfs(i);  
} 
//***************************************tarjan算法求强连通分量********************************//

void read(){
	freopen("data.txt","r",stdin);
	cout<<"start"<<endl; 
	int x,y,cur = 0;
	vector<pair<int,int> >E;
	vector<int>v;
	int same = 0; 
	while (scanf("%d%d",&x,&y)!=-1){
		//if (x == y && x < 1000)cout<<x<<endl;
		if (x == y)same++;
		E.push_back(make_pair(x,y));
		v.push_back(x);
		v.push_back(y);
	}
	cout<<"same:"<<same<<endl;
	sort(v.begin(),v.end());
	v.erase(unique(v.begin(), v.end()),v.end());
	tot = v.size();

	cout<<v[0]<<' '<<v[tot-1]<<endl;
	for (int i = 0; i < tot; ++i)
		pre_to_new_NUM[v[i]] = i;
	for (int i = 0; i < E.size(); ++i)
	{
		int xx = pre_to_new_NUM[E[i].first];
		int yy = pre_to_new_NUM[E[i].second];
		GG[xx].push_back(yy);
	}
	cout<<"start tarjan"<<endl;
	find_scc(tot);//tarjan
	cout<<"end tarjan"<<endl;	
	int maxx = 0;
	vector<Edge> edgeVec;
	vector<Node> nodeVec;
	edgemax = 0;
	for (int i = 0; i < tot; ++i){
		Node x;
		x.set_id(to_string(i));
		x.set_level(1);
		
	}
	graph.setNodesBundle(nodeVec);
	
	for (int i = 0; i < tot; ++i)
	for (int j = 0; j < GG[i].size(); ++j){
		int x = sccno[i] - 1;
		int y = sccno[GG[i][j]] - 1;
		if (x != y){
			G[x].push_back(y);
			Edge nm;
			nm.set_id(to_string(x)+"-"+to_string(y));
			nm.set_from(to_string(x));
			nm.set_to(to_string(y));
			edgeVec.push_back(nm);
		}
	}
	graph.setEdgesBundle(edgeVec);
	tot = scc_cnt;
	cout<<"tot="<<tot<<endl;
	
	
}
void toposort(){
	static int deg[maxn]={0};
	int cnt = 0;
	for (int i = 0; i < tot; ++i)
	for (int j = 0; j < G[i].size(); ++j){
		deg[G[i][j]]++;
		++cnt;
	} 
	cout<<cnt<<endl;
	cnt = 0;
	for (int i = 0; i < tot; ++i)
	{
		if (deg[i] != 0)cnt++;
	}
	printf("%d\n",cnt);
	queue<int>q;
	int sum = 0;
	for (int i = 0; i < tot; ++i)
	{
		if (deg[i] == 0){
			q.push(i);
			tln[i] = 1;
			++sum;
		}
	}
	int maxx = 0;
	cnt = 0;
	while (!q.empty()){
		++cnt;
		int x = q.front();
		q.pop();
		for (int i = 0; i < G[x].size(); ++i)
		{
			int y = G[x][i];
			--deg[y];
			tln[y] = max(tln[y], tln[x] + 1);
			maxx = max(tln[y], maxx);
			if (deg[y] == 0){
				q.push(y);
				sum += tln[y];
			}
		}
	}
	TF = (int)(log(maxx)/log(2)) + 1;
	cout<<cnt<<endl;
	cout<<"average:"<<sum*1.0/tot<<endl;
	cout<<"maxx:"<<maxx<<endl;
	cnt = 0;
	for (int i = 0; i < tot; ++i)
	{
		if (deg[i] > 0)cnt++;
	}
	printf("%d\n",cnt);
}

set<int> label_in[maxn*4], label_out[maxn*4];

int findfa(string s){
	int x = 0;
	//cout<<"##"<<endl;
	int i = 0;
	int len = s.length();
	while (i < len){
		if (s[i] >= '0' && s[i] <= '9'){
			x *= 10;
			x += s[i] - '0';
		}else break;
		++i;
	}
	//cout<<x<<endl;
	return x;	
} 

void construct(){	
	for (int i = 0; i < tot; ++i)
	{
		fa[i] = i;
	}
	auto start = graph.parse("match (a) return a"_graphsql); // return an iterator
	for (auto it=start; it!=graph.end(); ++it) //iterates over result set
	{
		Node a = it->getNode("a");
		a.set_level(tln[atio(a.id());
		graph.setNode(a);
	}
	for (int ii = 1; ii <= TF; ++ii){
	//a***********************************************************a
		auto start = graph.parse("match (a) return a"_graphsql); // return an iterator
		for (auto it=start; it!=graph.end(); ++it) //iterates over result set
		{
			Node a = it->getNode("a");
			int curl = a.level();
			if (curl % 2 == 0)continue;
		    set<string> outEdges = graph.getOutEdges(a.id());
		    vector<string> edge_remove;
		    bool flag = 0;
		    Node n;
		    for (string &e: outEdges)
		    {
		    	Node m = graph.getEdge(e).to();
		        if (m.level() > curl + 1)
		        {
		        	if (!flag){
		        		flag = 1;
		        		n = a;
						string name = a.id() + "#";
		        		n.set_id(name);
		        		n.set_level(a.level() + 1);
		        		graph.setNode(n);		        		
		        	}
		        	edge_remove.push_back(e);
					Edge nm;
					nm.set_id(a.id()+"-"+n.id());
					nm.set_from(a.id());
					nm.set_to(n.id());
					graph.setEdge(nm);
					
					nm.set_id(n.id()+"-"+m.id());
					nm.set_from(n.id());
					nm.set_to(m.id());
					graph.setEdge(nm);
		        }
		    }
		    for (int i = 0; i < edge_remove.size(); ++i)
		    	graph.removeEdge(edge_remove[i]);
        }
        //b**********************************************************************b
		auto start = graph.parse("match (a) return a"_graphsql); // return an iterator
		for (auto it=start; it!=graph.end(); ++it) //iterates over result set
		{
			Node a = it->getNode("a");
			int curl = a.level();
			if (curl % 2 == 0)continue;
		    set<string> inEdges = graph.getInEdges(a.id());
		    vector<string> edge_remove;
		    bool flag = 0;
		    Node n;
		    for (string &e: inEdges)
		    {
		    	Node m = graph.getEdge(e).from();
		    	int tmp = m.level();
		        if (tmp < curl - 1 && tmp % 2 == 0)
		        {
		        	if (!flag){
		        		flag = 1;
		        		n = a;
						string name = a.id() + "#";
		        		n.set_id(name);
		        		n.set_level(a.level() - 1);
		        		graph.setNode(n);		        		
		        	}
		        	edge_remove.push_back(e);
					Edge nm;
					nm.set_id(n.id() + "-" + a.id());
					nm.set_from(n.id());
					nm.set_to(a.id());
					graph.setEdge(nm);
					
					nm.set_id(m.id() + "-" + n.id());
					nm.set_from(m.id());
					nm.set_to(n.id());
					graph.setEdge(nm);
		        }
		    }
		    for (int i = 0; i < edge_remove.size(); ++i)
		    	graph.removeEdge(edge_remove[i]);
        }
    //c*****************************************************************c
		auto start = graph.parse("match (a) return a"_graphsql); // return an iterator
		for (auto it=start; it!=graph.end(); ++it) //iterates over result set
		{
			Node a = it->getNode("a");
			int curl = a.level();
			if (curl % 2 == 0)continue;
		    set<string> inEdges = graph.getInEdges(a.id());			
		    set<string> outEdges = graph.getOutEdges(a.id());
			set<string> inNode;
			set<string> outNode;
			int x = findfa(a.id());
		    for (string &e: outEdges)
		    {
		    	string to = graph.getEdge(e).to().id();
		    	outNode.insert(to);
				label_in[x].insert(findfa(to));		
			}
		    for (string &e: intEdges)
		    {
		    	string from = graph.getEdge(e).from.id();
		    	inNode.insert(from);			
		    	label_out[x].insert(findfa(from));
			}
			//��label 
			vector<Edge>edgeVec;
			for (string &to: outNode){
				for (string &from: inNode){
					Edge nm;
					nm.set_id(from+"-"+to);
					nm.set_from(from);
					nm.set_to(to);
					edgeVec.push_back(nm);					
				}
			}
			graph.setEdgesBundle(edgeVec);
			graph.eraseNode(a.id());						
		}    
	}	
}

//void organize_label(){}



int main()
{
	read();
	cout<<"Input completed."<<endl;
	toposort();
	construct();
	//organize_label();
	
	
	
	
	
	/*
	int x,y;
	vector<int> v;
	int cnt = 0, minn = 10000, maxx = 0;
	while (scanf("%d%d",&x,&y)!=-1){
	v.push_back(x);
	v.push_back(y);
		maxx = max(maxx,x);
		maxx = max(maxx,y);
		minn = min(minn,min(x,y));
		++cnt;
	}
	sort(v.begin(),v.end());
	v.erase(unique(v.begin(), v.end()),v.end());
	cout<<v.size()<<endl;
	cout<<maxx<<' '<<minn<<' '<<cnt<<endl;
	*/
	
	return 0;
}
