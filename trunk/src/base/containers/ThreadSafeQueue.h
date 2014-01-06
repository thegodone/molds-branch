//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   //
// Copyright (C) 2012-2014 Katushiko Nishimra                             // 
//                                                                        // 
// This file is part of MolDS.                                            // 
//                                                                        // 
// MolDS is free software: you can redistribute it and/or modify          // 
// it under the terms of the GNU General Public License as published by   // 
// the Free Software Foundation, either version 3 of the License, or      // 
// (at your option) any later version.                                    // 
//                                                                        // 
// MolDS is distributed in the hope that it will be useful,               // 
// but WITHOUT ANY WARRANTY; without even the implied warranty of         // 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          // 
// GNU General Public License for more details.                           // 
//                                                                        // 
// You should have received a copy of the GNU General Public License      // 
// along with MolDS.  If not, see <http://www.gnu.org/licenses/>.         // 
//************************************************************************//
#ifndef INCLUDED_THREADSAFEQUEQUE
#define INCLUDED_THREADSAFEQUEQUE
#include<queue>
#include<boost/shared_ptr.hpp>
#include<boost/thread.hpp>
#include<boost/thread/condition.hpp>
namespace MolDS_base_containers{

// This Queue class is thread-safe

template <typename T>
class ThreadSafeQueue
{
public:
   ThreadSafeQueue(){}
   ~ThreadSafeQueue(){}

   void Push(const T& data){
      boost::mutex::scoped_lock lk(this->stateGuard);
      this->stdQueue.push(data);
      this->stateChange.notify_all();
   } 

   T FrontPop(){
      boost::mutex::scoped_lock lk(this->stateGuard);
      if(this->stdQueue.empty()){
         std::stringstream ss;
         ss << "naitive queue has no member\n";
         MolDS_base::MolDSException ex(ss.str());
         int info = 0;
         ex.SetKeyValue<int>(MolDS_base::EmptyQueue, info);
         throw ex;
      }
      T ret = this->stdQueue.front();
      this->stdQueue.pop();
      this->stateChange.notify_all();
      return ret;
   }
     
   int Size(){
     boost::mutex::scoped_lock lk(this->stateGuard);
     return this->stdQueue.size();
   }
   
   bool Empty(){
     boost::mutex::scoped_lock lk(this->stateGuard);
     return this->stdQueue.empty();
   }
private:
   std::queue<T> stdQueue;
   boost::mutex  stateGuard;
   boost::condition_variable stateChange;
};

}
#endif
