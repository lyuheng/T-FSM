
#ifndef TASK_H_
#define TASK_H_

using namespace std;

template <class ContextT>
class Task
{
public:
	typedef ContextT ContextType; // used in comperOL.h
	ContextT context;

	// friend ibinstream &operator<<(ibinstream &m, const Task &t)
	// {
	// 	m << t.context;
	// 	return m;
	// }

	// friend obinstream &operator>>(obinstream &m, Task &t)
	// {
	// 	m >> t.context;
	// 	return m;
	// }

	// friend ifbinstream &operator<<(ifbinstream &m, const Task &t)
	// {
	// 	m << t.context;
	// 	return m;
	// }

	// friend ofbinstream &operator>>(ofbinstream &m, Task &t)
	// {
	// 	m >> t.context;
	// 	return m;
	// }
};

#endif /* TASK_H_ */
