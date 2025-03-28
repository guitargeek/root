// Author: Sergey Linev <s.linev@gsi.de>
// Date: 2017-10-16
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RWebWindow
#define ROOT7_RWebWindow

#include <ROOT/RWebDisplayHandle.hxx>

#include "ROOT/RConfig.hxx"

#include <memory>
#include <vector>
#include <string>
#include <queue>
#include <map>
#include <functional>
#include <mutex>
#include <thread>
#include <chrono>

class THttpCallArg;
class THttpServer;

namespace ROOT {

/// function signature for connect/disconnect call-backs
/// argument is connection id
using WebWindowConnectCallback_t = std::function<void(unsigned)>;

/// function signature for call-backs from the window clients
/// first argument is connection id, second is received data
using WebWindowDataCallback_t = std::function<void(unsigned, const std::string &)>;

/// function signature for waiting call-backs
/// Such callback used when calling thread need to waits for some special data,
/// but wants to run application event loop
/// As argument, spent time in second will be provided
/// Waiting will be performed until function returns non-zero value
using WebWindowWaitFunc_t = std::function<int(double)>;

class RFileDialog;
class RWebWindowsManager;
class RWebWindowWSHandler;

class RWebWindow {

   friend class RWebWindowsManager;
   friend class RWebWindowWSHandler;
   friend class RWebDisplayHandle;
   friend class RFileDialog;

private:
   using timestamp_t = std::chrono::time_point<std::chrono::system_clock>;

   struct QueueItem {
      int fChID{1};      ///<! channel
      bool fText{true};  ///<! is text data
      std::string fData; ///<! text or binary data
      QueueItem(int chid, bool txt, std::string &&data) : fChID(chid), fText(txt), fData(data) {}
   };

   struct WebConn {
      unsigned fConnId{0};                 ///<! connection id (unique inside the window)
      bool fHeadlessMode{false};           ///<! indicate if connection represent batch job
      bool fWasFirst{false};               ///<! indicate if this was first connection, will be reinjected also on first place
      std::string fKey;                    ///<! key value supplied to the window (when exists)
      int fKeyUsed{0};                     ///<! key value used to verify connection
      std::string fNewKey;                 ///<! new key if connection request reload
      std::unique_ptr<RWebDisplayHandle> fDisplayHandle;  ///<! handle assigned with started web display (when exists)
      std::shared_ptr<THttpCallArg> fHold; ///<! request used to hold headless browser
      timestamp_t fSendStamp;              ///<! last server operation, always used from window thread
      bool fActive{false};                 ///<! flag indicates if connection is active
      unsigned fWSId{0};                   ///<! websocket id
      int fReady{0};                       ///<! 0 - not ready, 1..9 - interim, 10 - done
      mutable std::mutex fMutex;           ///<! mutex must be used to protect all following data
      timestamp_t fRecvStamp;              ///<! last receive operation, protected with connection mutex
      int fRecvCount{0};                   ///<! number of received packets, should return back with next sending
      int fSendCredits{0};                 ///<! how many send operation can be performed without confirmation from other side
      int fClientCredits{0};               ///<! number of credits received from client
      bool fDoingSend{false};              ///<! true when performing send operation
      unsigned long fRecvSeq{0};           ///<! sequence id of last received packet
      unsigned long fSendSeq{1};           ///<! sequence id of last send packet
      std::queue<QueueItem> fQueue;        ///<! output queue
      std::map<int,std::shared_ptr<RWebWindow>> fEmbed; ///<! map of embed window for that connection, key value is channel id
      WebConn() = default;
      WebConn(unsigned connid) : fConnId(connid) {}
      WebConn(unsigned connid, unsigned wsid) : fConnId(connid), fActive(true), fWSId(wsid) {}
      WebConn(unsigned connid, bool headless_mode, const std::string &key)
         : fConnId(connid), fHeadlessMode(headless_mode), fKey(key)
      {
         ResetStamps();
      }
      ~WebConn();

      void ResetStamps() { fSendStamp = fRecvStamp = std::chrono::system_clock::now(); }

      void ResetData()
      {
         fActive = false;
         fWSId = 0;
         fReady = 0;
         fDoingSend = false;
         fSendCredits = 0;
         fClientCredits = 0;
         fRecvSeq = 0;
         fSendSeq = 1;
         while (!fQueue.empty())
            fQueue.pop();
      }
   };

   struct MasterConn {
      unsigned connid{0};
      int channel{-1};
      MasterConn(unsigned _connid, int _channel) : connid(_connid), channel(_channel) {}
   };

   enum EQueueEntryKind { kind_None, kind_Connect, kind_Data, kind_Disconnect };

   struct QueueEntry {
      unsigned fConnId{0};               ///<! connection id
      EQueueEntryKind fKind{kind_None};  ///<! kind of data
      std::string fData;                 ///<! data for given connection
      QueueEntry() = default;
      QueueEntry(unsigned connid, EQueueEntryKind kind, std::string &&data) : fConnId(connid), fKind(kind), fData(data) {}
   };

   using ConnectionsList_t = std::vector<std::shared_ptr<WebConn>>;

   std::shared_ptr<RWebWindowsManager> fMgr;        ///<! display manager
   std::shared_ptr<RWebWindow> fMaster;             ///<! master window where this window is embedded
   std::vector<MasterConn> fMasterConns;            ///<! master connections
   std::string fDefaultPage;                        ///<! HTML page (or file name) returned when window URL is opened
   std::string fPanelName;                          ///<! panel name which should be shown in the window
   unsigned fId{0};                                 ///<! unique identifier
   bool fUseServerThreads{false};                   ///<! indicates that server thread is using, no special window thread
   bool fUseProcessEvents{false};                   ///<! all window functionality will run through process events
   bool fProcessMT{false};                          ///<! if window event processing performed in dedicated thread
   bool fSendMT{false};                             ///<! true is special threads should be used for sending data
   bool fRequireAuthKey{true};                      ///<! defines if authentication key always required when connect to the widget
   std::shared_ptr<RWebWindowWSHandler> fWSHandler; ///<! specialize websocket handler for all incoming connections
   unsigned fConnCnt{0};                            ///<! counter of new connections to assign ids
   ConnectionsList_t fPendingConn;                  ///<! list of pending connection with pre-assigned keys
   ConnectionsList_t fConn;                         ///<! list of all accepted connections
   mutable std::mutex fConnMutex;                   ///<! mutex used to protect connection list
   unsigned fConnLimit{1};                          ///<! number of allowed active connections
   std::string fConnToken;                          ///<! value of "token" URL parameter which should be provided for connecting window
   bool fNativeOnlyConn{false};                     ///<! only native connection are allowed, created by Show() method
   bool fUseCurrentDir{false};                      ///<! if window can access local files via currentdir/ path of http server
   unsigned fMaxQueueLength{10};                    ///<! maximal number of queue entries
   WebWindowConnectCallback_t fConnCallback;        ///<! callback for connect event
   WebWindowDataCallback_t fDataCallback;           ///<! main callback when data over channel 1 is arrived
   WebWindowConnectCallback_t fDisconnCallback;     ///<! callback for disconnect event
   std::thread::id fCallbacksThrdId;                ///<! thread id where callbacks should be invoked
   bool fCallbacksThrdIdSet{false};                 ///<! flag indicating that thread id is assigned
   bool fHasWindowThrd{false};                      ///<! indicate if special window thread was started
   std::thread fWindowThrd;                         ///<! special thread for that window
   std::queue<QueueEntry> fInputQueue;              ///<! input queue for all callbacks
   std::mutex fInputQueueMutex;                     ///<! mutex to protect input queue
   unsigned fWidth{0}, fHeight{0};                  ///<! initial window width and height when displayed, zeros are ignored
   int fX{-1}, fY{-1};                              ///<! initial window position, -1 ignored
   float fOperationTmout{50.};                      ///<! timeout in seconds to perform synchronous operation, default 50s
   std::string fClientVersion;                      ///<! configured client version, used as prefix in scripts URL
   std::string fProtocolFileName;                   ///<! local file where communication protocol will be written
   int fProtocolCnt{-1};                            ///<! counter for protocol recording
   unsigned fProtocolConnId{0};                     ///<! connection id, which is used for writing protocol
   std::string fProtocolPrefix;                     ///<! prefix for created files names
   std::string fProtocol;                           ///<! protocol
   std::string fUserArgs;                           ///<! arbitrary JSON code, which is accessible via conn.getUserArgs() method
   std::shared_ptr<void> fClearOnClose;             ///<! entry which is cleared when last connection is closed
   static std::string gJSROOTsettings;              ///<! custom settings for JSROOT

   std::shared_ptr<RWebWindowWSHandler> CreateWSHandler(std::shared_ptr<RWebWindowsManager> mgr, unsigned id, double tmout);

   bool ProcessWS(THttpCallArg &arg);

   void CompleteWSSend(unsigned wsid);

   ConnectionsList_t GetWindowConnections(unsigned connid = 0, bool only_active = false) const;

   /// Find connection with specified websocket id
   std::shared_ptr<WebConn> FindConnection(unsigned wsid);

   void ClearConnection(std::shared_ptr<WebConn> &conn, bool provide_signal = false);

   std::shared_ptr<WebConn> RemoveConnection(unsigned wsid, bool provide_signal = false);

   bool _CanTrustIn(std::shared_ptr<WebConn> &conn, const std::string &key, const std::string &ntry, bool remote, bool test_first_time);

   std::string _MakeSendHeader(std::shared_ptr<WebConn> &conn, bool txt, const std::string &data, int chid);

   void ProvideQueueEntry(unsigned connid, EQueueEntryKind kind, std::string &&arg);

   void InvokeCallbacks(bool force = false);

   void SubmitData(unsigned connid, bool txt, std::string &&data, int chid = 1);

   bool CheckDataToSend(std::shared_ptr<WebConn> &conn);

   void CheckDataToSend(bool only_once = false);

   bool HasKey(const std::string &key, bool also_newkey = false) const;

   void RemoveKey(const std::string &key);

   std::string GenerateKey() const;

   void CheckPendingConnections();

   void CheckInactiveConnections();

   unsigned AddDisplayHandle(bool headless_mode, const std::string &key, std::unique_ptr<RWebDisplayHandle> &handle);

   unsigned AddEmbedWindow(std::shared_ptr<RWebWindow> window, unsigned connid, int channel);

   void RemoveEmbedWindow(unsigned connid, int channel);

   void AddMasterConnection(std::shared_ptr<RWebWindow> window, unsigned connid, int channel);

   std::vector<MasterConn> GetMasterConnections(unsigned connid = 0) const;

   void RemoveMasterConnection(unsigned connid = 0);

   bool ProcessBatchHolder(std::shared_ptr<THttpCallArg> &arg);

   std::string GetConnToken() const;

   unsigned MakeHeadless(bool create_new = false);

   unsigned FindHeadlessConnection();

   static std::function<bool(const std::shared_ptr<RWebWindow> &, unsigned, const std::string &)> gStartDialogFunc;

   static void SetStartDialogFunc(std::function<bool(const std::shared_ptr<RWebWindow> &, unsigned, const std::string &)>);

   static std::string HMAC(const std::string &key, const std::string &sessionKey, const char *msg, int msglen);

public:

   RWebWindow();

   ~RWebWindow();

   /// Returns ID for the window - unique inside window manager
   unsigned GetId() const { return fId; }

   /// Returns window manager
   std::shared_ptr<RWebWindowsManager> GetManager() const { return fMgr; }

   /// Set content of default window HTML page
   /// This page returns when URL address of the window will be requested
   /// Either HTML code or file name in the form "file:/home/user/data/file.htm"
   /// One also can using default locations like "file:rootui5sys/canv/canvas.html"
   void SetDefaultPage(const std::string &page) { fDefaultPage = page; }

   void SetPanelName(const std::string &name);

   /// Set window geometry. Will be applied if supported by used web display (like CEF or Chromium)
   void SetGeometry(unsigned width, unsigned height)
   {
      fWidth = width;
      fHeight = height;
   }

   /// Set window position. Will be applied if supported by used web display (like CEF or Chromium)
   void SetPosition(unsigned x, unsigned y)
   {
      fX = x;
      fY = y;
   }

   /////////////////////////////////////////////////////////////////////////
   /// returns configured window width (0 - default)
   /// actual window width can be different
   unsigned GetWidth() const { return fWidth; }

   /////////////////////////////////////////////////////////////////////////
   /// returns configured window height (0 - default)
   unsigned GetHeight() const { return fHeight; }

   /////////////////////////////////////////////////////////////////////////
   /// returns configured window X position (-1 - default)
   int GetX() const { return fX; }

   /////////////////////////////////////////////////////////////////////////
   /// returns configured window Y position (-1 - default)
   int GetY() const { return fY; }

   void SetConnLimit(unsigned lmt = 0);

   unsigned GetConnLimit() const;

   void SetConnToken(const std::string &token = "");

   /////////////////////////////////////////////////////////////////////////
   /// configures maximal queue length of data which can be held by window
   void SetMaxQueueLength(unsigned len = 10) { fMaxQueueLength = len; }

   /////////////////////////////////////////////////////////////////////////
   /// Return maximal queue length of data which can be held by window
   unsigned GetMaxQueueLength() const { return fMaxQueueLength; }

   /////////////////////////////////////////////////////////////////////////
   /// configures that only native (own-created) connections are allowed
   void SetNativeOnlyConn(bool on = true) { fNativeOnlyConn = on; }

   /////////////////////////////////////////////////////////////////////////
   /// returns true if only native (own-created) connections are allowed
   bool IsNativeOnlyConn() const { return fNativeOnlyConn; }

   /////////////////////////////////////////////////////////////////////////
   /// Configure if authentication key in connection string is required
   void SetRequireAuthKey(bool on) { fRequireAuthKey = on; }

   /////////////////////////////////////////////////////////////////////////
   /// returns true if authentication string is required
   bool IsRequireAuthKey() const { return fRequireAuthKey; }

   /////////////////////////////////////////////////////////////////////////
   /// Configure if window can access local files via currentdir/ path of http server
   void SetUseCurrentDir(bool on = true) { fUseCurrentDir = on; }

   /////////////////////////////////////////////////////////////////////////
   /// returns true if window can access local files via currentdir/ path of http server
   bool IsUseCurrentDir() const { return fUseCurrentDir; }

   void SetClientVersion(const std::string &vers);

   std::string GetClientVersion() const;

   void SetUserArgs(const std::string &args);

   std::string GetUserArgs() const;

   int NumConnections(bool with_pending = false) const;

   unsigned GetConnectionId(int num = 0) const;

   std::vector<unsigned> GetConnections(unsigned excludeid = 0) const;

   bool HasConnection(unsigned connid = 0, bool only_active = true) const;

   void CloseConnections();

   void CloseConnection(unsigned connid);

   /// Returns timeout for synchronous WebWindow operations
   float GetOperationTmout() const { return fOperationTmout; }

   /// Set timeout for synchronous WebWindow operations
   void SetOperationTmout(float tm = 50.) { fOperationTmout = tm; }

   std::string GetUrl(bool remote = true);

   THttpServer *GetServer();

   void Sync();

   void Run(double tm = 0.);

   unsigned Show(const RWebDisplayArgs &args = "");

   unsigned GetDisplayConnection() const;

   /// Returns true when window was shown at least once
   bool IsShown() const { return GetDisplayConnection() != 0; }

   bool CanSend(unsigned connid, bool direct = true) const;

   int GetSendQueueLength(unsigned connid) const;

   void Send(unsigned connid, const std::string &data);

   void SendBinary(unsigned connid, const void *data, std::size_t len);

   void SendBinary(unsigned connid, std::string &&data);

   void RecordData(const std::string &fname = "protocol.json", const std::string &fprefix = "");

   std::string GetAddr() const;

   _R__DEPRECATED_LATER("Use GetUrl() to get valid connection URL") std::string GetRelativeAddr(const std::shared_ptr<RWebWindow> &win) const;

   _R__DEPRECATED_LATER("Use GetAddr() to get valid connection URL") std::string GetRelativeAddr(const RWebWindow &win) const;

   void SetCallBacks(WebWindowConnectCallback_t conn, WebWindowDataCallback_t data, WebWindowConnectCallback_t disconn = nullptr);

   void Reset();

   void SetConnectCallBack(WebWindowConnectCallback_t func);

   void SetDataCallBack(WebWindowDataCallback_t func);

   void SetDisconnectCallBack(WebWindowConnectCallback_t func);

   void SetClearOnClose(const std::shared_ptr<void> &handle = nullptr);

   void AssignThreadId();

   void UseServerThreads();

   int WaitFor(WebWindowWaitFunc_t check);

   int WaitForTimed(WebWindowWaitFunc_t check);

   int WaitForTimed(WebWindowWaitFunc_t check, double duration);

   void StartThread();

   void StopThread();

   void TerminateROOT();

   static std::shared_ptr<RWebWindow> Create();

   static unsigned ShowWindow(std::shared_ptr<RWebWindow> window, const RWebDisplayArgs &args = "");

   static bool IsFileDialogMessage(const std::string &msg);

   static bool EmbedFileDialog(const std::shared_ptr<RWebWindow> &window, unsigned connid, const std::string &args);

   static void SetJSROOTSettings(const std::string &json);
};

} // namespace ROOT

#endif
