(window.webpackJsonp = window.webpackJsonp || []).push([[0], {
  10: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__(15),
      __webpack_require__(30),
      __webpack_require__(24),
      __webpack_require__(7),
      __webpack_require__(21),
      __webpack_require__(5),
      __webpack_require__(26);
      var react = __webpack_require__(0)
        , react_default = __webpack_require__.n(react)
        , prop_types = __webpack_require__(1)
        , prop_types_default = __webpack_require__.n(prop_types)
        , memoize = __webpack_require__(464)
        , omit = __webpack_require__(924)
        , forOwn = __webpack_require__(914)
        , PropTypes = __webpack_require__(57)
        , es = __webpack_require__(127)
        , reduce = (__webpack_require__(35),
      __webpack_require__(927))
        , pick = __webpack_require__(319)
        , isEqual = __webpack_require__(922)
        , createMSAStore = __webpack_require__(230)
        , actions = __webpack_require__(8)
        , positionReducers = __webpack_require__(47)
        , requestAnimation = __webpack_require__(153)
        , conservation_worker = __webpack_require__(229)
        , conservation_worker_default = __webpack_require__.n(conservation_worker);
      function _extends() {
          return (_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      function _objectWithoutProperties(source, excluded) {
          if (null == source)
              return {};
          var key, i, target = function _objectWithoutPropertiesLoose(source, excluded) {
              if (null == source)
                  return {};
              var key, i, target = {}, sourceKeys = Object.keys(source);
              for (i = 0; i < sourceKeys.length; i++)
                  key = sourceKeys[i],
                  0 <= excluded.indexOf(key) || (target[key] = source[key]);
              return target
          }(source, excluded);
          if (Object.getOwnPropertySymbols) {
              var sourceSymbolKeys = Object.getOwnPropertySymbols(source);
              for (i = 0; i < sourceSymbolKeys.length; i++)
                  key = sourceSymbolKeys[i],
                  0 <= excluded.indexOf(key) || Object.prototype.propertyIsEnumerable.call(source, key) && (target[key] = source[key])
          }
          return target
      }
      let worker = null;
      const setUpWorker = (store,sequences,sampleSize=null,element)=>{
          worker = new conservation_worker_default.a,
          worker.postMessage({
              sequences: sequences,
              sampleSize: sampleSize
          }),
          worker.onmessage = e=>{
              store.dispatch(actions.b.updateConservation(e.data)),
              element && element.current && element.current.el && element.current.el.current && element.current.el.current.dispatchEvent(new CustomEvent("conservationProgress",{
                  bubbles: !0,
                  detail: e.data
              })),
              1 === e.data.progress && (console.log("completed conservation analisys"),
              store.dispatch(actions.b.updateSequences(sequences)))
          }
      }
        , reduxActions = {
          sequences: "updateSequences"
      };
      Object.keys(PropTypes.b).forEach(key=>{
          !(key in reduxActions) && PropTypes.b[key] && (reduxActions[key] = "updateProp")
      }
      );
      const attributesToStore = Object.keys(reduxActions)
        , mapToActionKeys = obj=>Object(reduce.a)(obj, (acc,v)=>(acc[v.key] = v,
      acc), {})
        , mainStoreActionKeys = mapToActionKeys(actions.b)
        , positionStoreActionKeys = mapToActionKeys(positionReducers.a);
      var propsToRedux_ref = react_default.a.createElement("div", null, " Error initializing the MSAViewer. ");
      var propsToRedux = WrappedComponent=>{
          class PropsToReduxComponent extends react.Component {
              constructor(props) {
                  super(props);
                  const storeProps = Object(pick.a)(props, attributesToStore) || {};
                  this.el = react_default.a.createRef(),
                  this.msaStore = props.msaStore,
                  void 0 !== storeProps.sequences ? (this.msaStore = Object(createMSAStore.a)(storeProps),
                  storeProps.calculateConservation && this.msaStore && setUpWorker(this.msaStore, storeProps.sequences, storeProps.sampleSizeConservation, this.el)) : console.warn("Check your MSA properties", storeProps)
              }
              componentDidMount() {
                  void 0 !== this.props.position && this.updatePosition(this.props.position)
              }
              componentDidUpdate(oldProps) {
                  const newProps = this.props;
                  for (const prop in Object(pick.a)(newProps, attributesToStore))
                      if (!Object(isEqual.a)(oldProps[prop], newProps[prop]))
                          if ("position" === prop)
                              this.updatePosition(newProps[prop]);
                          else if (prop in reduxActions) {
                              let action;
                              switch ("calculateConservation" === prop && (newProps[prop] ? setUpWorker(this.msaStore, this.props.sequences, this.props.sampleSizeConservation, this.el) : worker.terminate()),
                              "sampleSizeConservation" === prop && (worker && worker.terminate(),
                              setUpWorker(this.msaStore, this.props.sequences, this.props.sampleSizeConservation, this.el)),
                              reduxActions[prop]) {
                              case "updateProp":
                                  action = actions.b[reduxActions[prop]](prop, newProps[prop]);
                                  break;
                              default:
                                  action = actions.b[reduxActions[prop]](newProps[prop])
                              }
                              this.msaStore.dispatch(action)
                          } else
                              console.error(prop, " is unknown.")
              }
              dispatch(action) {
                  if (action.type in mainStoreActionKeys)
                      this.msaStore.dispatch(action);
                  else {
                      if (!(action.type in positionStoreActionKeys))
                          throw new Error("Invalid action",action);
                      this.el.current.positionStore.dispatch(action)
                  }
              }
              getColorMap() {
                  const colorScheme = this.props.colorScheme;
                  let map = {};
                  try {
                      map = this.msaStore.getState().props.colorScheme.scheme.map
                  } catch (_unused) {}
                  return {
                      name: colorScheme,
                      map: map
                  }
              }
              render() {
                  const _omit = Object(omit.a)(this.props, attributesToStore)
                    , msaStore = _omit.msaStore
                    , props = _objectWithoutProperties(_omit, ["msaStore"]);
                  return void 0 === this.msaStore ? propsToRedux_ref : react_default.a.createElement(WrappedComponent, _extends({
                      ref: this.el,
                      msaStore: msaStore || this.msaStore
                  }, props))
              }
          }
          return PropsToReduxComponent.displayName = "PropsToReduxComponent",
          Object(forOwn.a)(actions.b, (v,k)=>{
              PropsToReduxComponent.prototype[k] = function(payload) {
                  this.msaStore.dispatch(v(payload))
              }
          }
          ),
          Object(forOwn.a)(positionReducers.a, (v,k)=>{
              PropsToReduxComponent.prototype[k] = function(payload) {
                  Object(requestAnimation.a)(this, ()=>{
                      this.el.current.positionStore.dispatch(v(payload))
                  }
                  )
              }
          }
          ),
          PropsToReduxComponent
      }
        , redux = __webpack_require__(138);
      const same = "FORWARD_SAME_PROP_NAME";
      function forwardPropsMapper(props, propsToForward) {
          const forward = {};
          let remainingProps = props;
          return Object(forOwn.a)(propsToForward, (v,k)=>{
              const result = function util_forwardProps(props, propsSelector) {
                  const forward = {}
                    , other = {};
                  return Object(forOwn.a)(props, (v,k)=>{
                      if (k in propsSelector) {
                          const forwardedName = propsSelector[k];
                          forward[forwardedName === same ? k : forwardedName] = v
                      } else
                          other[k] = v
                  }
                  ),
                  {
                      forward: forward,
                      other: other
                  }
              }(remainingProps, v);
              forward[k] = result.forward,
              remainingProps = result.other
          }
          ),
          {
              forward: forward,
              other: remainingProps
          }
      }
      var Labels = __webpack_require__(254)
        , OverviewBar = __webpack_require__(255)
        , PositionBar = __webpack_require__(131)
        , SequenceViewer = __webpack_require__(97)
        , some = __webpack_require__(930)
        , shallowEqual = __webpack_require__(93);
      class PureBaseLayout_PureBaseLayout extends react.Component {
          shouldComponentUpdate(nextProps) {
              return !(Object(shallowEqual.a)(this.props, nextProps) && Object(some.a)(this.props.forwardedPropsKeys.map(k=>Object(shallowEqual.a)(this.props[k], nextProps[k]))))
          }
      }
      PureBaseLayout_PureBaseLayout.propTypes = {
          forwardPropKeys: prop_types_default.a.array
      };
      var layouts_PureBaseLayout = PureBaseLayout_PureBaseLayout
        , connect = __webpack_require__(16);
      function ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function _defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      function basic_extends() {
          return (basic_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      class basic_MSABasicLayout extends layouts_PureBaseLayout {
          render() {
              const labelsPadding = this.props.tileHeight;
              return react_default.a.createElement("div", {
                  style: {
                      display: "flex"
                  }
              }, react_default.a.createElement(Labels.a, basic_extends({
                  style: {
                      paddingTop: labelsPadding + 50
                  }
              }, this.props.labelsProps)), react_default.a.createElement("div", null, react_default.a.createElement(OverviewBar.a, basic_extends({
                  height: 50
              }, this.props.overviewBarProps)), react_default.a.createElement(PositionBar.a, this.props.positionBarProps), react_default.a.createElement(SequenceViewer.a, this.props.sequenceViewerProps), react_default.a.createElement("div", {
                  style: {
                      height: 10
                  }
              })))
          }
      }
      basic_MSABasicLayout.displayName = "MSABasicLayout",
      basic_MSABasicLayout.propTypes = function _objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                  _defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }({}, layouts_PureBaseLayout.propTypes);
      basic_MSABasicLayout.__docgenInfo = {
          description: "",
          methods: [],
          displayName: "MSABasicLayout",
          composes: ["./PureBaseLayout"]
      };
      var basic = Object(connect.a)(state=>({
          tileHeight: state.props.tileHeight
      }))(basic_MSABasicLayout);
      function inverse_ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function inverse_defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      function inverse_extends() {
          return (inverse_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/layouts/basic.js"] = {
          name: "MSABasicLayout",
          docgenInfo: basic_MSABasicLayout.__docgenInfo,
          path: "src/components/layouts/basic.js"
      });
      class inverse_MSAInverseLayout extends layouts_PureBaseLayout {
          render() {
              const labelsPadding = this.props.tileHeight;
              return react_default.a.createElement("div", {
                  style: {
                      display: "flex"
                  }
              }, react_default.a.createElement("div", null, react_default.a.createElement(OverviewBar.a, inverse_extends({
                  height: 50
              }, this.props.overviewBarProps)), react_default.a.createElement(PositionBar.a, this.props.positionBarProps), react_default.a.createElement(SequenceViewer.a, this.props.sequenceViewerProps)), react_default.a.createElement(Labels.a, inverse_extends({
                  style: {
                      paddingTop: labelsPadding + 50
                  }
              }, this.props.labelsProps)))
          }
      }
      inverse_MSAInverseLayout.displayName = "MSAInverseLayout",
      inverse_MSAInverseLayout.propTypes = function inverse_objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? inverse_ownKeys(Object(source), !0).forEach((function(key) {
                  inverse_defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : inverse_ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }({}, layouts_PureBaseLayout.propTypes);
      inverse_MSAInverseLayout.__docgenInfo = {
          description: "",
          methods: [],
          displayName: "MSAInverseLayout",
          composes: ["./PureBaseLayout"]
      };
      var inverse = Object(connect.a)(state=>({
          tileHeight: state.props.tileHeight
      }))(inverse_MSAInverseLayout);
      function full_extends() {
          return (full_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/layouts/inverse.js"] = {
          name: "MSAInverseLayout",
          docgenInfo: inverse_MSAInverseLayout.__docgenInfo,
          path: "src/components/layouts/inverse.js"
      });
      var full_ref = react_default.a.createElement(SequenceViewer.a, null)
        , _ref2 = react_default.a.createElement("br", null)
        , _ref3 = react_default.a.createElement("br", null)
        , _ref4 = react_default.a.createElement(OverviewBar.a, null);
      class full_MSAFullLayout extends layouts_PureBaseLayout {
          render() {
              return react_default.a.createElement("div", {
                  style: {
                      display: "flex"
                  }
              }, react_default.a.createElement(Labels.a, this.props.labelsProps), react_default.a.createElement("div", null, full_ref, react_default.a.createElement(PositionBar.a, this.props.positionBarProps), _ref2, react_default.a.createElement(OverviewBar.a, this.props.overviewBarProps), _ref3, react_default.a.createElement(PositionBar.a, this.props.positionBarProps), _ref4, react_default.a.createElement(OverviewBar.a, full_extends({}, this.props.overviewBarProps, {
                  method: "information-content"
              }))))
          }
      }
      full_MSAFullLayout.displayName = "MSAFullLayout",
      full_MSAFullLayout.__docgenInfo = {
          description: "",
          methods: [],
          displayName: "MSAFullLayout"
      };
      var full = full_MSAFullLayout;
      function compact_ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function compact_defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/layouts/full.js"] = {
          name: "MSAFullLayout",
          docgenInfo: full_MSAFullLayout.__docgenInfo,
          path: "src/components/layouts/full.js"
      });
      class compact_MSACompactLayout extends layouts_PureBaseLayout {
          render() {
              return react_default.a.createElement("div", null, react_default.a.createElement(PositionBar.a, this.props.positionBarProps), react_default.a.createElement(SequenceViewer.a, this.props.sequenceViewerProps))
          }
      }
      compact_MSACompactLayout.displayName = "MSACompactLayout",
      compact_MSACompactLayout.propTypes = function compact_objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? compact_ownKeys(Object(source), !0).forEach((function(key) {
                  compact_defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : compact_ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }({}, layouts_PureBaseLayout.propTypes),
      compact_MSACompactLayout.__docgenInfo = {
          description: "",
          methods: [],
          displayName: "MSACompactLayout",
          composes: ["./PureBaseLayout"]
      };
      var compact = compact_MSACompactLayout;
      function funky_ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function funky_defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      function funky_extends() {
          return (funky_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/layouts/compact.js"] = {
          name: "MSACompactLayout",
          docgenInfo: compact_MSACompactLayout.__docgenInfo,
          path: "src/components/layouts/compact.js"
      });
      var funky_ref = react_default.a.createElement("br", null);
      class funky_MSAFunkyLayout extends layouts_PureBaseLayout {
          render() {
              const labelsStyle = {
                  paddingTop: this.props.tileHeight
              };
              return react_default.a.createElement("div", {
                  style: {
                      display: "flex"
                  }
              }, react_default.a.createElement(Labels.a, {
                  style: labelsStyle
              }), react_default.a.createElement("div", null, react_default.a.createElement(PositionBar.a, this.props.positionBarProps), react_default.a.createElement(SequenceViewer.a, this.props.sequenceViewerProps), react_default.a.createElement(PositionBar.a, this.props.positionBarProps), react_default.a.createElement(OverviewBar.a, this.props.overviewBarProps), funky_ref, react_default.a.createElement(PositionBar.a, this.props.positionBarProps)), react_default.a.createElement(Labels.a, funky_extends({
                  style: labelsStyle
              }, this.props.labelsProps)))
          }
      }
      funky_MSAFunkyLayout.displayName = "MSAFunkyLayout",
      funky_MSAFunkyLayout.propTypes = function funky_objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? funky_ownKeys(Object(source), !0).forEach((function(key) {
                  funky_defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : funky_ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }({}, layouts_PureBaseLayout.propTypes);
      funky_MSAFunkyLayout.__docgenInfo = {
          description: "",
          methods: [],
          displayName: "MSAFunkyLayout",
          composes: ["./PureBaseLayout"]
      };
      var funky = Object(connect.a)(state=>({
          tileHeight: state.props.tileHeight
      }))(funky_MSAFunkyLayout);
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/layouts/funky.js"] = {
          name: "MSAFunkyLayout",
          docgenInfo: funky_MSAFunkyLayout.__docgenInfo,
          path: "src/components/layouts/funky.js"
      });
      var isEmpty = __webpack_require__(923)
        , withPositionStore = (__webpack_require__(178),
      __webpack_require__(37));
      function Coordinates_extends() {
          return (Coordinates_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      function Coordinates_ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function Coordinates_objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? Coordinates_ownKeys(Object(source), !0).forEach((function(key) {
                  Coordinates_defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : Coordinates_ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }
      function Coordinates_defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      function Coordinates_objectWithoutProperties(source, excluded) {
          if (null == source)
              return {};
          var key, i, target = function Coordinates_objectWithoutPropertiesLoose(source, excluded) {
              if (null == source)
                  return {};
              var key, i, target = {}, sourceKeys = Object.keys(source);
              for (i = 0; i < sourceKeys.length; i++)
                  key = sourceKeys[i],
                  0 <= excluded.indexOf(key) || (target[key] = source[key]);
              return target
          }(source, excluded);
          if (Object.getOwnPropertySymbols) {
              var sourceSymbolKeys = Object.getOwnPropertySymbols(source);
              for (i = 0; i < sourceSymbolKeys.length; i++)
                  key = sourceSymbolKeys[i],
                  0 <= excluded.indexOf(key) || Object.prototype.propertyIsEnumerable.call(source, key) && (target[key] = source[key])
          }
          return target
      }
      const Coordinate = _ref=>{
          let coordinateComponent = _ref.coordinateComponent
            , _ref$coordinateAttrib = _ref.coordinateAttributes
            , coordinateAttributes = void 0 === _ref$coordinateAttrib ? {} : _ref$coordinateAttrib
            , _ref$coordinateStyle = _ref.coordinateStyle
            , coordinateStyle = void 0 === _ref$coordinateStyle ? {} : _ref$coordinateStyle
            , otherProps = Coordinates_objectWithoutProperties(_ref, ["coordinateComponent", "coordinateAttributes", "coordinateStyle"]);
          if (coordinateComponent)
              return react_default.a.createElement(coordinateComponent, otherProps);
          {
              const start = otherProps.start
                , end = otherProps.end
                , tileHeight = otherProps.tileHeight
                , attributes = Coordinates_objectSpread({}, otherProps, {}, coordinateAttributes);
              return react_default.a.createElement("div", Coordinates_extends({
                  style: Coordinates_objectSpread({
                      height: tileHeight
                  }, coordinateStyle)
              }, attributes), start, "-", end)
          }
      }
      ;
      class Coordinates_Coordinates extends react.PureComponent {
          constructor(props) {
              super(props),
              this.el = react_default.a.createRef()
          }
          shouldRerender() {
              return !0
          }
          render() {
              const _this$props = this.props
                , height = _this$props.height
                , width = _this$props.width
                , tileHeight = _this$props.tileHeight
                , tileWidth = _this$props.tileWidth
                , position = _this$props.position
                , coordinateComponent = _this$props.coordinateComponent
                , coordinateAttributes = _this$props.coordinateAttributes
                , coordinateStyle = _this$props.coordinateStyle
                , sequences = _this$props.sequences;
              this.props.position.lastCurrentViewSequence = position.currentViewSequence,
              this.props.position.lastStartYTile = 0;
              const start = Math.round(position.xPos / tileWidth) + 1
                , end = start + Math.round(width / tileWidth) - 1
                , firstSequenceInView = position.currentViewSequence
                , nSequencesInView = Math.ceil(height / tileHeight)
                , spacer = react_default.a.createElement("div", {
                  style: {
                      height: firstSequenceInView * tileHeight
                  }
              });
              return react_default.a.createElement("div", null, react_default.a.createElement("div", {
                  style: {
                      height: height,
                      overflow: "hidden",
                      position: "relative",
                      whiteSpace: "nowrap"
                  },
                  ref: this.el
              }, spacer, sequences.slice(position.currentViewSequence, firstSequenceInView + nSequencesInView + 1).map((sequence,index)=>react_default.a.createElement(Coordinate, {
                  key: index,
                  start: start,
                  end: end,
                  tileHeight: tileHeight,
                  coordinateComponent: coordinateComponent,
                  coordinateStyle: coordinateStyle,
                  coordinateAttributes: coordinateAttributes,
                  sequence: sequence
              }))))
          }
      }
      Coordinates_Coordinates.displayName = "Coordinates",
      Coordinates_Coordinates.propTypes = {
          width: prop_types_default.a.number.isRequired,
          height: prop_types_default.a.number.isRequired,
          tileWidth: prop_types_default.a.number.isRequired,
          tileHeight: prop_types_default.a.number.isRequired,
          position: prop_types_default.a.shape({
              currentViewSequence: prop_types_default.a.number.isRequired,
              currentViewSequencePosition: prop_types_default.a.number.isRequired,
              lastCurrentViewSequence: prop_types_default.a.number,
              lastStartYTile: prop_types_default.a.number,
              xPos: prop_types_default.a.number.isRequired,
              xPosOffset: prop_types_default.a.number.isRequired,
              yPos: prop_types_default.a.number.isRequired,
              yPosOffset: prop_types_default.a.number.isRequired
          }).isRequired,
          coordinateComponent: prop_types_default.a.oneOfType([prop_types_default.a.object, prop_types_default.a.func]),
          coordinateStyle: prop_types_default.a.object,
          coordinateAttributes: prop_types_default.a.object
      },
      Coordinates_Coordinates.defaultProps = {
          coordinateComponent: void 0,
          coordinateStyle: {},
          coordinateAttributes: {}
      };
      Coordinates_Coordinates.__docgenInfo = {
          description: "Displays the coordinates of the current scrolled position.",
          methods: [{
              name: "shouldRerender",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }],
          displayName: "Coordinates",
          props: {
              coordinateComponent: {
                  defaultValue: {
                      value: "undefined",
                      computed: !0
                  },
                  type: {
                      name: "union",
                      value: [{
                          name: "object"
                      }, {
                          name: "func"
                      }]
                  },
                  required: !1,
                  description: "Component to create labels from."
              },
              coordinateStyle: {
                  defaultValue: {
                      value: "{}",
                      computed: !1
                  },
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Inline styles to apply to each label."
              },
              coordinateAttributes: {
                  defaultValue: {
                      value: "{}",
                      computed: !1
                  },
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Attributes to apply to each label."
              },
              width: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: "Width of the sequence viewer (in pixels), e.g. `500`."
              },
              height: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: "Height of the sequence viewer (in pixels), e.g. `500`."
              },
              tileWidth: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: "Width of the main tiles (in pixels), e.g. `20`"
              },
              tileHeight: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: "Height of the main tiles (in pixels), e.g. `20`"
              },
              position: {
                  type: {
                      name: "shape",
                      value: {
                          currentViewSequence: {
                              name: "number",
                              required: !0
                          },
                          currentViewSequencePosition: {
                              name: "number",
                              required: !0
                          },
                          lastCurrentViewSequence: {
                              name: "number",
                              required: !1
                          },
                          lastStartYTile: {
                              name: "number",
                              required: !1
                          },
                          xPos: {
                              name: "number",
                              required: !0
                          },
                          xPosOffset: {
                              name: "number",
                              required: !0
                          },
                          yPos: {
                              name: "number",
                              required: !0
                          },
                          yPosOffset: {
                              name: "number",
                              required: !0
                          }
                      }
                  },
                  required: !0,
                  description: ""
              }
          }
      };
      var yBars_Coordinates = Object(connect.a)(state=>({
          width: state.props.width,
          height: state.props.height,
          tileHeight: state.props.tileHeight,
          tileWidth: state.props.tileWidth,
          sequences: state.sequences.raw
      }))(Object(withPositionStore.a)(Coordinates_Coordinates, {
          withX: !0,
          withY: !0
      }));
      function nightingale_ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function nightingale_defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/yBars/Coordinates.js"] = {
          name: "Coordinates",
          docgenInfo: Coordinates_Coordinates.__docgenInfo,
          path: "src/components/yBars/Coordinates.js"
      });
      class nightingale_NightingaleLayout extends layouts_PureBaseLayout {
          render() {
              const _this$props = this.props
                , leftCoordinatesProps = _this$props.leftCoordinatesProps
                , rightCoordinatesProps = _this$props.rightCoordinatesProps
                , labelsProps = _this$props.labelsProps;
              return react_default.a.createElement("div", {
                  style: {
                      display: "flex"
                  }
              }, !Object(isEmpty.a)(labelsProps) && react_default.a.createElement(Labels.a, labelsProps), !Object(isEmpty.a)(leftCoordinatesProps) && react_default.a.createElement(yBars_Coordinates, leftCoordinatesProps), react_default.a.createElement(SequenceViewer.a, this.props.sequenceViewerProps), !Object(isEmpty.a)(rightCoordinatesProps) && react_default.a.createElement(yBars_Coordinates, rightCoordinatesProps))
          }
      }
      nightingale_NightingaleLayout.displayName = "NightingaleLayout",
      nightingale_NightingaleLayout.propTypes = function nightingale_objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? nightingale_ownKeys(Object(source), !0).forEach((function(key) {
                  nightingale_defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : nightingale_ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }({}, layouts_PureBaseLayout.propTypes),
      nightingale_NightingaleLayout.__docgenInfo = {
          description: "",
          methods: [],
          displayName: "NightingaleLayout",
          composes: ["./PureBaseLayout"]
      };
      var nightingale = nightingale_NightingaleLayout;
      function MSALayouter_extends() {
          return (MSALayouter_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      function MSALayouter_objectWithoutProperties(source, excluded) {
          if (null == source)
              return {};
          var key, i, target = function MSALayouter_objectWithoutPropertiesLoose(source, excluded) {
              if (null == source)
                  return {};
              var key, i, target = {}, sourceKeys = Object.keys(source);
              for (i = 0; i < sourceKeys.length; i++)
                  key = sourceKeys[i],
                  0 <= excluded.indexOf(key) || (target[key] = source[key]);
              return target
          }(source, excluded);
          if (Object.getOwnPropertySymbols) {
              var sourceSymbolKeys = Object.getOwnPropertySymbols(source);
              for (i = 0; i < sourceSymbolKeys.length; i++)
                  key = sourceSymbolKeys[i],
                  0 <= excluded.indexOf(key) || Object.prototype.propertyIsEnumerable.call(source, key) && (target[key] = source[key])
          }
          return target
      }
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/layouts/nightingale.js"] = {
          name: "NightingaleLayout",
          docgenInfo: nightingale_NightingaleLayout.__docgenInfo,
          path: "src/components/layouts/nightingale.js"
      });
      const layouts = {
          basic: basic,
          inverse: inverse,
          full: full,
          compact: compact,
          funky: funky,
          nightingale: nightingale,
          default: basic
      };
      class MSALayouter_MSALayouter extends react.PureComponent {
          constructor(props) {
              super(props),
              this.el = react_default.a.createRef(),
              this.forwardedPropsKeys = Object.keys(this.constructor.propsToForward)
          }
          dispatchEvent(event) {
              this.el.current.dispatchEvent(event)
          }
          render() {
              const _this$props = this.props
                , layout = _this$props.layout
                , otherProps = MSALayouter_objectWithoutProperties(_this$props, ["layout"]);
              if (layout in layouts) {
                  const Layout = layouts[layout]
                    , _forwardPropsMapper = forwardPropsMapper(otherProps, this.constructor.propsToForward)
                    , forwardProps = _forwardPropsMapper.forward
                    , propsOnDiv = _forwardPropsMapper.other;
                  return react_default.a.createElement("div", MSALayouter_extends({
                      el: this.ref
                  }, propsOnDiv), react_default.a.createElement(Layout, MSALayouter_extends({}, forwardProps, {
                      forwardedPropsKeys: this.forwardedPropsKeys
                  })))
              }
              console.error("$this.props.layout} is invalid. Please use one of ".concat(Object.keys(layouts)))
          }
      }
      MSALayouter_MSALayouter.displayName = "MSALayouter",
      function MSALayouter_defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }(MSALayouter_MSALayouter, "propsToForward", {
          sequenceViewerProps: {
              showModBar: same,
              onResidueMouseEnter: same,
              onResidueMouseLeave: same,
              onResidueClick: same,
              onResidueDoubleClick: same,
              highlight: same,
              features: same,
              onFeatureClick: same,
              onHighlightMouseEnter: same,
              onHighlightMouseLeave: same,
              sequenceBorder: "border",
              sequenceBorderColor: "borderColor",
              sequenceBorderWidth: "borderWidth",
              sequenceTextColor: "textColor",
              sequenceTextFont: "textFont",
              sequenceOverflow: "overflow",
              sequenceOverflowX: "overflowX",
              sequenceOverflowy: "overflowY",
              sequenceScrollBarPositionX: "scrollBarPositionX",
              sequenceScrollBarPositionY: "scrollBarPositionY",
              sequenceDisableDragging: same,
              overlayConservation: same,
              onPositionUpdate: same
          },
          leftCoordinatesProps: {
              leftCoordinate: same,
              leftCoordinateComponent: "coordinateComponent",
              leftCoordinateStyle: "coordinateStyle",
              leftCoordinateAttributes: "coordinateAttributes"
          },
          rightCoordinatesProps: {
              rightCoordinate: same,
              rightCoordinateComponent: "coordinateComponent",
              rightCoordinateStyle: "coordinateStyle",
              rightCoordinateAttributes: "coordinateAttributes"
          },
          labelsProps: {
              labelComponent: same,
              labelStyle: same,
              labelAttributes: same
          },
          positionBarProps: {
              markerSteps: same,
              markerStartIndex: "startIndex",
              markerComponent: same,
              markerStyle: same,
              markerAttributes: same
          },
          overviewBarProps: {
              barMethod: "method",
              barFillColor: "fillColor",
              barComponent: same,
              barStyle: same,
              barAttributes: same
          }
      }),
      MSALayouter_MSALayouter.defaultProps = {
          layout: "default"
      },
      MSALayouter_MSALayouter.propTyes = {
          layout: prop_types_default.a.oneOf(Object.keys(layouts))
      },
      MSALayouter_MSALayouter.__docgenInfo = {
          description: "Pick the selected layout and forwards all properties to it.",
          methods: [{
              name: "dispatchEvent",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "event",
                  type: null
              }],
              returns: null
          }],
          displayName: "MSALayouter",
          props: {
              layout: {
                  defaultValue: {
                      value: '"default"',
                      computed: !1
                  },
                  required: !1
              }
          }
      };
      var layouts_MSALayouter = MSALayouter_MSALayouter;
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/layouts/MSALayouter.js"] = {
          name: "MSALayouter",
          docgenInfo: MSALayouter_MSALayouter.__docgenInfo,
          path: "src/components/layouts/MSALayouter.js"
      });
      var shallowSelect = __webpack_require__(76);
      function MSAViewer_extends() {
          return (MSAViewer_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      function MSAViewer_objectWithoutProperties(source, excluded) {
          if (null == source)
              return {};
          var key, i, target = function MSAViewer_objectWithoutPropertiesLoose(source, excluded) {
              if (null == source)
                  return {};
              var key, i, target = {}, sourceKeys = Object.keys(source);
              for (i = 0; i < sourceKeys.length; i++)
                  key = sourceKeys[i],
                  0 <= excluded.indexOf(key) || (target[key] = source[key]);
              return target
          }(source, excluded);
          if (Object.getOwnPropertySymbols) {
              var sourceSymbolKeys = Object.getOwnPropertySymbols(source);
              for (i = 0; i < sourceSymbolKeys.length; i++)
                  key = sourceSymbolKeys[i],
                  0 <= excluded.indexOf(key) || Object.prototype.propertyIsEnumerable.call(source, key) && (target[key] = source[key])
          }
          return target
      }
      function MSAViewer_ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function MSAViewer_objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? MSAViewer_ownKeys(Object(source), !0).forEach((function(key) {
                  MSAViewer_defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : MSAViewer_ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }
      function MSAViewer_defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      const defaultEvents = {
          onResidueClick: "residueClick"
      };
      class MSAViewer_MSAViewerComponent extends react.Component {
          constructor(props) {
              super(props),
              this.el = react_default.a.createRef(),
              this._setupStores(),
              this.createDomHandler = Object(memoize.a)(this.createDomHandler.bind(this)),
              this.forwardProps = Object(shallowSelect.a)(p=>Object(omit.a)(p, ["msaStore"]), this.forwardProps.bind(this))
          }
          _setupStores() {
              this.positionStore = Object(redux.b)(positionReducers.c),
              this.positionStore.dispatch(positionReducers.a.updateMainStore(this.props.msaStore.getState())),
              this.msaStoreUnsubscribe = this.props.msaStore.subscribe(()=>{
                  this.positionStore.dispatch(positionReducers.a.updateMainStore(this.props.msaStore.getState()))
              }
              )
          }
          componentWillUnmount() {
              this.msaStoreUnsubscribe()
          }
          getChildContext() {
              return {
                  positionMSAStore: this.positionStore
              }
          }
          createDomHandler(domEventName) {
              return e=>{
                  const event = new CustomEvent(domEventName,{
                      detail: e,
                      bubbles: !0
                  });
                  this.el.current.dispatchEvent(event)
              }
          }
          forwardProps(props) {
              const options = MSAViewer_objectSpread({}, props);
              return Object(forOwn.a)(defaultEvents, (forwardedName,currentName)=>{
                  currentName in defaultEvents && !(currentName in options) && (options[currentName] = this.createDomHandler(forwardedName))
              }
              ),
              options
          }
          render() {
              const _this$props = this.props
                , children = _this$props.children
                , msaStore = _this$props.msaStore
                , otherProps = MSAViewer_objectWithoutProperties(_this$props, ["children", "msaStore"]);
              return children ? react_default.a.createElement(es.a, {
                  store: msaStore
              }, react_default.a.createElement("div", MSAViewer_extends({}, otherProps, {
                  ref: this.el
              }), children)) : react_default.a.createElement(es.a, {
                  store: msaStore
              }, react_default.a.createElement("div", {
                  ref: this.el
              }, react_default.a.createElement(layouts_MSALayouter, this.forwardProps(otherProps))))
          }
      }
      MSAViewer_MSAViewerComponent.displayName = "MSAViewerComponent",
      MSAViewer_MSAViewerComponent.childContextTypes = {
          positionMSAStore: prop_types_default.a.object
      },
      MSAViewer_MSAViewerComponent.propTypes = {
          sequences: prop_types_default.a.arrayOf(PropTypes.d).isRequired,
          width: prop_types_default.a.number,
          height: prop_types_default.a.number,
          tileWidth: prop_types_default.a.number,
          tileHeight: prop_types_default.a.number,
          position: PropTypes.c,
          colorScheme: PropTypes.a,
          onResidueMouseEnter: prop_types_default.a.func,
          onResidueMouseLeave: prop_types_default.a.func,
          onResidueClick: prop_types_default.a.func,
          onResidueDoubleClick: prop_types_default.a.func,
          onFeatureClick: prop_types_default.a.func,
          highlight: prop_types_default.a.oneOfType([prop_types_default.a.object, prop_types_default.a.array]),
          features: prop_types_default.a.arrayOf(prop_types_default.a.object),
          layout: prop_types_default.a.oneOf(["basic", "default", "inverse", "full", "compact", "funky", "nightingale"]),
          sequenceBorder: prop_types_default.a.bool,
          sequenceBorderColor: prop_types_default.a.string,
          sequenceBorderWidth: prop_types_default.a.number,
          sequenceTextColor: prop_types_default.a.string,
          sequenceTextFont: prop_types_default.a.string,
          sequenceOverflow: prop_types_default.a.oneOf(["hidden", "auto", "scroll"]),
          sequenceOverflowX: prop_types_default.a.oneOf(["hidden", "auto", "scroll", "initial"]),
          sequenceOverflowY: prop_types_default.a.oneOf(["hidden", "auto", "scroll", "initial"]),
          sequenceScrollBarPositionX: prop_types_default.a.oneOf(["top", "bottom"]),
          sequenceScrollBarPositionY: prop_types_default.a.oneOf(["left", "right"]),
          labelComponent: prop_types_default.a.oneOfType([prop_types_default.a.object, prop_types_default.a.func]),
          labelStyle: prop_types_default.a.object,
          labelAttributes: prop_types_default.a.object,
          markerSteps: prop_types_default.a.number,
          markerStartIndex: prop_types_default.a.number,
          markerComponent: prop_types_default.a.oneOfType([prop_types_default.a.object, prop_types_default.a.func]),
          markerStyle: prop_types_default.a.object,
          markerAttributes: prop_types_default.a.object,
          barMethod: prop_types_default.a.oneOf(["information-content", "conservation"]),
          barFillColor: prop_types_default.a.string,
          barStyle: prop_types_default.a.object,
          barAttributes: prop_types_default.a.object,
          msaStore: prop_types_default.a.object
      };
      const MSAViewer = propsToRedux(MSAViewer_MSAViewerComponent);
      MSAViewer.defaultProps = PropTypes.e,
      MSAViewer.propTypes = MSAViewer_objectSpread({}, MSAViewer_MSAViewerComponent.propTypes),
      MSAViewer_MSAViewerComponent.propTypes = Object(omit.a)(MSAViewer_MSAViewerComponent.propTypes, [...Object.keys(PropTypes.b), "sequences", "position"]);
      __webpack_exports__.a = MSAViewer;
      MSAViewer_MSAViewerComponent.__docgenInfo = {
          description: "A general-purpose layout for the MSA components\n\nWhen children are passed it acts as a Context Provider for the msaStore,\notherwise it provides a default layout and forwards it props the respective\ncomponents.",
          methods: [{
              name: "_setupStores",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "createDomHandler",
              docblock: "Creates a listener which triggers `domEventName`",
              modifiers: [],
              params: [{
                  name: "domEventName"
              }],
              returns: null,
              description: "Creates a listener which triggers `domEventName`"
          }, {
              name: "forwardProps",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "props",
                  type: null
              }],
              returns: null
          }],
          displayName: "MSAViewerComponent",
          props: {
              sequences: {
                  type: {
                      name: "arrayOf",
                      value: {
                          name: "custom",
                          raw: "SequencePropType"
                      }
                  },
                  required: !0,
                  description: 'Sequence data.\n`sequences` expects an array of individual sequences.\n\n`sequence`: Raw sequence, e.g. `MEEPQSDPSIEP` (required)\n`name`: name of the sequence, e.g. `Sequence X`\n\nExample:\n\n```js\nconst sequences = [\n  {\n    name: "seq.1",\n    sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED",\n  },\n  {\n    name: "seq.2",\n    sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n];\n```'
              },
              width: {
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Width of the sequence viewer (in pixels), e.g. `500`."
              },
              height: {
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Height of the sequence viewer (in pixels), e.g. `500`."
              },
              tileWidth: {
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Width of the main tiles (in pixels), e.g. `20`"
              },
              tileHeight: {
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Height of the main tiles (in pixels), e.g. `20`"
              },
              position: {
                  type: {
                      name: "custom",
                      raw: "PositionPropType"
                  },
                  required: !1,
                  description: "Current x and y position of the viewpoint\nin the main sequence viewer (in pixels).\nThis specifies the position of the top-left corner\nof the viewpoint within the entire alignment,\ne.g. `{xPos: 20, yPos: 5}`."
              },
              colorScheme: {
                  type: {
                      name: "custom",
                      raw: "ColorSchemePropType"
                  },
                  required: !1,
                  description: "Colorscheme to use. Currently the follow colorschemes are supported:\n`aliphatic`, `aromatic`, `buried`, `buried_index`, `charged`, `cinema`, `clustal2`,\n`clustal`, `helix`, `helix_propensity`, `hydro`, `lesk`, `mae`, `negative`,\n`nucleotide`, `polar`, `positive`, `purine`, `purine_pyrimidine`, `serine_threonine`,\n`strand`, `strand_propensity`, `taylor`, `turn`, `turn_propensity`, `zappo`, `conservation`\n\nSee [msa-colorschemes](https://github.com/wilzbach/msa-colorschemes) for details."
              },
              onResidueMouseEnter: {
                  type: {
                      name: "func"
                  },
                  required: !1,
                  description: "Callback fired when the mouse pointer is entering a residue."
              },
              onResidueMouseLeave: {
                  type: {
                      name: "func"
                  },
                  required: !1,
                  description: "Callback fired when the mouse pointer is leaving a residue."
              },
              onResidueClick: {
                  type: {
                      name: "func"
                  },
                  required: !1,
                  description: "Callback fired when the mouse pointer clicked a residue."
              },
              onResidueDoubleClick: {
                  type: {
                      name: "func"
                  },
                  required: !1,
                  description: "Callback fired when the mouse pointer clicked a residue."
              },
              onFeatureClick: {
                  type: {
                      name: "func"
                  },
                  required: !1,
                  description: "Callback fired when the mouse pointer clicked a feature."
              },
              highlight: {
                  type: {
                      name: "union",
                      value: [{
                          name: "object"
                      }, {
                          name: "array"
                      }]
                  },
                  required: !1,
                  description: "Displays a highlight"
              },
              features: {
                  type: {
                      name: "arrayOf",
                      value: {
                          name: "object"
                      }
                  },
                  required: !1,
                  description: "An array of features which can be clicked"
              },
              layout: {
                  type: {
                      name: "enum",
                      value: [{
                          value: '"basic"',
                          computed: !1
                      }, {
                          value: '"default"',
                          computed: !1
                      }, {
                          value: '"inverse"',
                          computed: !1
                      }, {
                          value: '"full"',
                          computed: !1
                      }, {
                          value: '"compact"',
                          computed: !1
                      }, {
                          value: '"funky"',
                          computed: !1
                      }, {
                          value: '"nightingale"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: "Predefined layout scheme to use (only used when no child elements are provided).\nAvailable layouts: `basic`, `inverse`, `full`, `compact`, `funky`, `nightingale`"
              },
              sequenceBorder: {
                  type: {
                      name: "bool"
                  },
                  required: !1,
                  description: "Whether to draw a border."
              },
              sequenceBorderColor: {
                  type: {
                      name: "string"
                  },
                  required: !1,
                  description: "Color of the border. Name, hex or RGB value."
              },
              sequenceBorderWidth: {
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Width of the border."
              },
              sequenceTextColor: {
                  type: {
                      name: "string"
                  },
                  required: !1,
                  description: "Color of the text residue letters (name, hex or RGB value)"
              },
              sequenceTextFont: {
                  type: {
                      name: "string"
                  },
                  required: !1,
                  description: "Font to use when drawing the individual residues."
              },
              sequenceOverflow: {
                  type: {
                      name: "enum",
                      value: [{
                          value: '"hidden"',
                          computed: !1
                      }, {
                          value: '"auto"',
                          computed: !1
                      }, {
                          value: '"scroll"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: "What should happen if content overflows."
              },
              sequenceOverflowX: {
                  type: {
                      name: "enum",
                      value: [{
                          value: '"hidden"',
                          computed: !1
                      }, {
                          value: '"auto"',
                          computed: !1
                      }, {
                          value: '"scroll"',
                          computed: !1
                      }, {
                          value: '"initial"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: 'What should happen if x-axis content overflows (overwrites "overflow")'
              },
              sequenceOverflowY: {
                  type: {
                      name: "enum",
                      value: [{
                          value: '"hidden"',
                          computed: !1
                      }, {
                          value: '"auto"',
                          computed: !1
                      }, {
                          value: '"scroll"',
                          computed: !1
                      }, {
                          value: '"initial"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: 'What should happen if y-axis content overflows (overwrites "overflow")'
              },
              sequenceScrollBarPositionX: {
                  type: {
                      name: "enum",
                      value: [{
                          value: '"top"',
                          computed: !1
                      }, {
                          value: '"bottom"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: 'X Position of the scroll bar ("top or "bottom")'
              },
              sequenceScrollBarPositionY: {
                  type: {
                      name: "enum",
                      value: [{
                          value: '"left"',
                          computed: !1
                      }, {
                          value: '"right"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: 'Y Position of the scroll bar ("left" or "right")'
              },
              labelComponent: {
                  type: {
                      name: "union",
                      value: [{
                          name: "object"
                      }, {
                          name: "func"
                      }]
                  },
                  required: !1,
                  description: "Component to create labels from."
              },
              labelStyle: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Inline styles to apply to each label."
              },
              labelAttributes: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Attributes to apply to each label."
              },
              markerSteps: {
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "At which steps the position labels should appear, e.g. `2` for (1, 3, 5)"
              },
              markerStartIndex: {
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "At which number the PositionBar marker should start counting.\nTypical values are: `1` (1-based indexing) and `0` (0-based indexing)."
              },
              markerComponent: {
                  type: {
                      name: "union",
                      value: [{
                          name: "object"
                      }, {
                          name: "func"
                      }]
                  },
                  required: !1,
                  description: "Component to create markers from."
              },
              markerStyle: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Inline styles to apply to each marker."
              },
              markerAttributes: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Attributes to apply to each marker."
              },
              barMethod: {
                  type: {
                      name: "enum",
                      value: [{
                          value: '"information-content"',
                          computed: !1
                      }, {
                          value: '"conservation"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: "Method to use for the OverviewBar:\n - `information-content`: Information entropy after Shannon of a column (scaled)\n - `conservation`: Conservation of a column (scaled)"
              },
              barFillColor: {
                  type: {
                      name: "string"
                  },
                  required: !1,
                  description: "Fill color of the OverviewBar, e.g. `#999999`"
              },
              barStyle: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Inline styles to apply to each bar."
              },
              barAttributes: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Attributes to apply to each bar."
              },
              msaStore: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "A custom msaStore (created with `createMSAStore`).\nUseful for custom interaction with other components"
              }
          },
          childContext: {
              positionMSAStore: {
                  type: {
                      name: "object"
                  },
                  required: !1
              }
          }
      },
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/MSAViewer.js"] = {
          name: "MSAViewerComponent",
          docgenInfo: MSAViewer_MSAViewerComponent.__docgenInfo,
          path: "src/components/MSAViewer.js"
      })
  },
  112: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.d(__webpack_exports__, "a", (function() {
          return ColorScheme
      }
      )),
      __webpack_require__.d(__webpack_exports__, "b", (function() {
          return isColorScheme
      }
      ));
      const schemesMgr = new (__webpack_require__(187).a);
      class ColorScheme {
          constructor(colorScheme) {
              this.scheme = schemesMgr.getScheme(colorScheme)
          }
          updateConservation(conservation) {
              schemesMgr.conservation = conservation
          }
          getColor(element, position) {
              return this.scheme.getColor(element, position, schemesMgr.conservation)
          }
      }
      function isColorScheme(obj) {
          return obj && "function" == typeof obj.getColor
      }
  },
  117: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      const assert = __webpack_require__(891);
      __webpack_exports__.a = assert
  },
  131: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__(15),
      __webpack_require__(30),
      __webpack_require__(7),
      __webpack_require__(21),
      __webpack_require__(5);
      var react__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(0)
        , react__WEBPACK_IMPORTED_MODULE_5___default = __webpack_require__.n(react__WEBPACK_IMPORTED_MODULE_5__)
        , prop_types__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(1)
        , prop_types__WEBPACK_IMPORTED_MODULE_6___default = __webpack_require__.n(prop_types__WEBPACK_IMPORTED_MODULE_6__)
        , lodash_es__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(915)
        , lodash_es__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(319)
        , _store_connect__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(16)
        , _utils_shallowSelect__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(76)
        , _utils_autobind__WEBPACK_IMPORTED_MODULE_11__ = __webpack_require__(77)
        , _xBar__WEBPACK_IMPORTED_MODULE_12__ = __webpack_require__(249);
      function _extends() {
          return (_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      function ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function _objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                  _defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }
      function _defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      function _objectWithoutProperties(source, excluded) {
          if (null == source)
              return {};
          var key, i, target = function _objectWithoutPropertiesLoose(source, excluded) {
              if (null == source)
                  return {};
              var key, i, target = {}, sourceKeys = Object.keys(source);
              for (i = 0; i < sourceKeys.length; i++)
                  key = sourceKeys[i],
                  0 <= excluded.indexOf(key) || (target[key] = source[key]);
              return target
          }(source, excluded);
          if (Object.getOwnPropertySymbols) {
              var sourceSymbolKeys = Object.getOwnPropertySymbols(source);
              for (i = 0; i < sourceSymbolKeys.length; i++)
                  key = sourceSymbolKeys[i],
                  0 <= excluded.indexOf(key) || Object.prototype.propertyIsEnumerable.call(source, key) && (target[key] = source[key])
          }
          return target
      }
      class HTMLPositionBarComponent extends react__WEBPACK_IMPORTED_MODULE_5__.PureComponent {
          constructor(props) {
              super(props),
              this.cache = function() {}
              ,
              Object(_utils_autobind__WEBPACK_IMPORTED_MODULE_11__.a)(this, "createMarker"),
              this.marker = Object(_utils_shallowSelect__WEBPACK_IMPORTED_MODULE_10__.a)(Object(lodash_es__WEBPACK_IMPORTED_MODULE_7__.a)(lodash_es__WEBPACK_IMPORTED_MODULE_8__.a, this.constructor.markerAttributes), this.createMarker)
          }
          createMarker(props) {
              return this.cache = function() {}
              ,
              function createMarker({markerSteps: markerSteps, startIndex: startIndex, tileWidth: tileWidth, font: font, markerComponent: markerComponent, markerStyle: markerStyle, markerAttributes: markerAttributes}) {
                  class Marker extends react__WEBPACK_IMPORTED_MODULE_5__.PureComponent {
                      render() {
                          const _this$props = this.props
                            , index = _this$props.index
                            , otherProps = _objectWithoutProperties(_this$props, ["index"]);
                          if (markerComponent)
                              return react__WEBPACK_IMPORTED_MODULE_5___default.a.createElement(markerComponent, {
                                  index: index
                              });
                          {
                              let name;
                              otherProps.style = _objectSpread({
                                  width: tileWidth,
                                  display: "inline-block",
                                  textAlign: "center"
                              }, markerStyle),
                              name = 0 == index % markerSteps ? index + 0 + startIndex : ".";
                              const attributes = _objectSpread({}, otherProps, {}, markerAttributes);
                              return react__WEBPACK_IMPORTED_MODULE_5___default.a.createElement("div", attributes, name)
                          }
                      }
                  }
                  return Marker.displayName = "Marker",
                  Marker
              }(props)
          }
          render() {
              const _this$props2 = this.props
                , cacheElements = _this$props2.cacheElements
                , otherProps = (_this$props2.markerSteps,
              _this$props2.startIndex,
              _this$props2.dispatch,
              _this$props2.markerComponent,
              _this$props2.markerStyle,
              _objectWithoutProperties(_this$props2, ["cacheElements", "markerSteps", "startIndex", "dispatch", "markerComponent", "markerStyle"]));
              return react__WEBPACK_IMPORTED_MODULE_5___default.a.createElement(_xBar__WEBPACK_IMPORTED_MODULE_12__.a, _extends({
                  tileComponent: this.marker(this.props),
                  cacheElements: cacheElements,
                  componentCache: this.cache
              }, otherProps))
          }
      }
      HTMLPositionBarComponent.displayName = "HTMLPositionBarComponent",
      _defineProperty(HTMLPositionBarComponent, "markerAttributes", ["markerSteps", "startIndex", "tileWidth", "markerComponent", "markerStyle", "markerAttributes"]),
      HTMLPositionBarComponent.defaultProps = {
          style: {
              font: "12px Arial"
          },
          height: 15,
          markerSteps: 2,
          startIndex: 1,
          cacheElements: 10,
          markerStyle: {}
      },
      HTMLPositionBarComponent.propTypes = {
          font: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.string,
          height: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.number,
          markerSteps: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.number,
          startIndex: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.number,
          markerComponent: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.oneOfType([prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.object, prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.func]),
          style: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.object,
          markerStyle: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.object,
          markerAttributes: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.object
      };
      HTMLPositionBarComponent.__docgenInfo = {
          description: "Displays the sequence names with an arbitrary Marker component",
          methods: [{
              name: "createMarker",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "props",
                  type: null
              }],
              returns: null
          }],
          displayName: "HTMLPositionBarComponent",
          props: {
              style: {
                  defaultValue: {
                      value: '{\n  font: "12px Arial",\n}',
                      computed: !1
                  },
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Inline styles to apply to the PositionBar component"
              },
              height: {
                  defaultValue: {
                      value: "15",
                      computed: !1
                  },
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Height of the PositionBar (in pixels), e.g. `100`"
              },
              markerSteps: {
                  defaultValue: {
                      value: "2",
                      computed: !1
                  },
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "At which steps the position labels should appear, e.g. `2` for (1, 3, 5)"
              },
              startIndex: {
                  defaultValue: {
                      value: "1",
                      computed: !1
                  },
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "At which number the PositionBar marker should start counting.\nTypical values are: `1` (1-based indexing) and `0` (0-based indexing)."
              },
              cacheElements: {
                  defaultValue: {
                      value: "10",
                      computed: !1
                  },
                  required: !1
              },
              markerStyle: {
                  defaultValue: {
                      value: "{}",
                      computed: !1
                  },
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Inline styles to apply to each marker."
              },
              font: {
                  type: {
                      name: "string"
                  },
                  required: !1,
                  description: "Font of the sequence labels, e.g. `20px Arial`"
              },
              markerComponent: {
                  type: {
                      name: "union",
                      value: [{
                          name: "object"
                      }, {
                          name: "func"
                      }]
                  },
                  required: !1,
                  description: "Component to create markers from."
              },
              markerAttributes: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Attributes to apply to each marker."
              }
          }
      },
      __webpack_exports__.a = Object(_store_connect__WEBPACK_IMPORTED_MODULE_9__.a)(state=>({
          sequences: state.sequences.raw,
          maxLength: state.sequences.maxLength,
          width: state.props.width,
          tileWidth: state.props.tileWidth,
          nrXTiles: state.sequenceStats.nrXTiles
      }))(HTMLPositionBarComponent),
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/xBars/PositionBar.js"] = {
          name: "HTMLPositionBarComponent",
          docgenInfo: HTMLPositionBarComponent.__docgenInfo,
          path: "src/components/xBars/PositionBar.js"
      })
  },
  153: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_exports__.a = function requestAnimation(instance, callback) {
          void 0 === instance.nextFrame && (instance.nextFrame = window.requestAnimationFrame(function() {
              callback(),
              this.nextFrame = void 0
          }
          .bind(instance)))
      }
  },
  158: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_exports__.a = !1
  },
  159: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      var react = __webpack_require__(0)
        , react_default = __webpack_require__.n(react)
        , prop_types = __webpack_require__(1)
        , prop_types_default = __webpack_require__.n(prop_types);
      var base = class DrawingBase {
          constructor(el) {
              this.el = el,
              this.state = {}
          }
          updateEl(el) {
              this.el = el
          }
          font(fontName) {
              this.state.font = fontName
          }
          globalAlpha(globalAlpha) {
              this.state.globalAlpha = globalAlpha
          }
          startDrawingFrame() {
              this.clear()
          }
          endDrawingFrame() {}
          save() {}
          restore() {}
      }
      ;
      var cache = class CanvasCharCache {
          constructor() {
              this.cache = {},
              this.cacheHeight = 0,
              this.cacheWidth = 0
          }
          getFontTile(letter, width, height, font) {
              return width === this.cacheWidth && height === this.cacheHeight && font === this.font || (this.updateDimensions(width, height),
              this.font = font),
              void 0 === this.cache[letter] && this.createTile(letter, width, height),
              this.cache[letter]
          }
          createTile(letter, width, height, font) {
              const canvas = this.cache[letter] = document.createElement("canvas");
              return canvas.width = width,
              canvas.height = height,
              this.ctx = canvas.getContext("2d"),
              this.ctx.font = this.font + "px mono",
              this.ctx.textBaseline = "middle",
              this.ctx.textAlign = "center",
              this.ctx.fillText(letter, width / 2, height / 2, width, font)
          }
          updateDimensions(width, height) {
              this.invalidate(),
              this.cacheWidth = width,
              this.cacheHeight = height
          }
          invalidate() {
              this.cache = {}
          }
      }
      ;
      var canvas = class canvas_Canvas extends base {
          constructor(el) {
              super(el),
              this.ctx = el.getContext("2d"),
              this.cache = new cache
          }
          clear() {
              this.ctx.clearRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height)
          }
          fillRect(x, y, width, height) {
              this.ctx.fillRect(x, y, width, height)
          }
          fillText(text, x, y, width, height) {
              return this.ctx.drawImage(this.cache.getFontTile(text, width, height, this.ctx.font), x, y, width, height)
          }
          font(fontName) {
              this.ctx.font = fontName
          }
          fillStyle(fillStyle) {
              this.ctx.fillStyle = fillStyle
          }
          globalAlpha(globalAlpha) {
              this.ctx.globalAlpha = globalAlpha
          }
          save() {
              this.ctx.save()
          }
          restore() {
              this.ctx.restore()
          }
      }
      ;
      class CanvasComponent_CanvasComponent extends react.PureComponent {
          constructor(props) {
              super(props),
              this.canvas = react_default.a.createRef()
          }
          componentDidMount() {
              this.ctx = new canvas(this.canvas.current),
              this.draw()
          }
          componentDidUpdate() {
              this._draw()
          }
          _draw() {
              this.ctx && (this.ctx.startDrawingFrame(),
              this.ctx.save(),
              this.draw(),
              this.ctx.restore(),
              this.ctx.endDrawingFrame())
          }
          draw() {
              console.error("Implement me.")
          }
          render() {
              return react_default.a.createElement("div", {
                  style: this.props.style
              }, react_default.a.createElement("canvas", {
                  ref: this.canvas,
                  width: this.props.width,
                  height: this.props.height
              }))
          }
      }
      CanvasComponent_CanvasComponent.displayName = "CanvasComponent",
      function _defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }(CanvasComponent_CanvasComponent, "defaultProps", {
          engine: "canvas"
      }),
      CanvasComponent_CanvasComponent.propTypes = {
          width: prop_types_default.a.number.isRequired,
          height: prop_types_default.a.number.isRequired,
          style: prop_types_default.a.object,
          engine: prop_types_default.a.oneOf(["canvas", "webgl"])
      },
      CanvasComponent_CanvasComponent.__docgenInfo = {
          description: "Constructs a drawable canvas (e.g. HTML Canvas or WebGL) and provides it as\na reference.\n\nOn every redraw, this.draw() gets called.",
          methods: [{
              name: "_draw",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "draw",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }],
          displayName: "CanvasComponent",
          props: {
              engine: {
                  defaultValue: {
                      value: '"canvas"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"canvas"',
                          computed: !1
                      }, {
                          value: '"webgl"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: "Rendering engine: `canvas` or `webgl` (experimental)."
              },
              width: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: "Width of the component (in pixels), e.g. `100`"
              },
              height: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: "Width of the component (in pixels), e.g. `100`"
              },
              style: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Custom style configuration."
              }
          }
      };
      __webpack_exports__.a = CanvasComponent_CanvasComponent;
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/Canvas/CanvasComponent.js"] = {
          name: "CanvasComponent",
          docgenInfo: CanvasComponent_CanvasComponent.__docgenInfo,
          path: "src/components/Canvas/CanvasComponent.js"
      })
  },
  16: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      var react_redux__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(127);
      __webpack_exports__.a = function msaConnect(mapStateToProps, mapDispatchToProps, mergeProps, options={}) {
          return Object(react_redux__WEBPACK_IMPORTED_MODULE_0__.b)(mapStateToProps, mapDispatchToProps, mergeProps, options)
      }
  },
  187: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.d(__webpack_exports__, "c", (function() {
          return staticSchemes
      }
      )),
      __webpack_require__.d(__webpack_exports__, "b", (function() {
          return dynSchemes
      }
      ));
      __webpack_require__(35);
      const StaticSchemeClass = function(map) {
          this.defaultColor = "#ffffff",
          this.type = "static",
          this.map = map,
          this.getColor = function(letter) {
              return void 0 !== this.map[letter] ? this.map[letter] : this.defaultColor
          }
      }
        , DynSchemeClass = function(fun, opt) {
          this.type = "dyn",
          this.opt = opt,
          void 0 !== fun.init ? (fun.init.call(this),
          this.getColor = fun.run,
          this.reset = fun.init,
          this.map = fun.map) : this.getColor = fun
      };
      var buried = {
          A: "#00a35c",
          R: "#00fc03",
          N: "#00eb14",
          D: "#00eb14",
          C: "#0000ff",
          Q: "#00f10e",
          E: "#00f10e",
          G: "#009d62",
          H: "#00d52a",
          I: "#0054ab",
          L: "#007b84",
          K: "#00ff00",
          M: "#009768",
          F: "#008778",
          P: "#00e01f",
          S: "#00d52a",
          T: "#00db24",
          W: "#00a857",
          Y: "#00e619",
          V: "#005fa0",
          B: "#00eb14",
          X: "#00b649",
          Z: "#00f10e"
      }
        , helix = {
          A: "#e718e7",
          R: "#6f906f",
          N: "#1be41b",
          D: "#778877",
          C: "#23dc23",
          Q: "#926d92",
          E: "#ff00ff",
          G: "#00ff00",
          H: "#758a75",
          I: "#8a758a",
          L: "#ae51ae",
          K: "#a05fa0",
          M: "#ef10ef",
          F: "#986798",
          P: "#00ff00",
          S: "#36c936",
          T: "#47b847",
          W: "#8a758a",
          Y: "#21de21",
          V: "#857a85",
          B: "#49b649",
          X: "#758a75",
          Z: "#c936c9"
      }
        , purine = {
          A: " #FF83FA",
          C: " #40E0D0",
          G: " #FF83FA",
          R: " #FF83FA",
          T: " #40E0D0",
          U: " #40E0D0",
          Y: " #40E0D0"
      }
        , strand = {
          A: "#5858a7",
          R: "#6b6b94",
          N: "#64649b",
          D: "#2121de",
          C: "#9d9d62",
          Q: "#8c8c73",
          E: "#0000ff",
          G: "#4949b6",
          H: "#60609f",
          I: "#ecec13",
          L: "#b2b24d",
          K: "#4747b8",
          M: "#82827d",
          F: "#c2c23d",
          P: "#2323dc",
          S: "#4949b6",
          T: "#9d9d62",
          W: "#c0c03f",
          Y: "#d3d32c",
          V: "#ffff00",
          B: "#4343bc",
          X: "#797986",
          Z: "#4747b8"
      }
        , turn = {
          A: "#2cd3d3",
          R: "#708f8f",
          N: "#ff0000",
          D: "#e81717",
          C: "#a85757",
          Q: "#3fc0c0",
          E: "#778888",
          G: "#ff0000",
          H: "#708f8f",
          I: "#00ffff",
          L: "#1ce3e3",
          K: "#7e8181",
          M: "#1ee1e1",
          F: "#1ee1e1",
          P: "#f60909",
          S: "#e11e1e",
          T: "#738c8c",
          W: "#738c8c",
          Y: "#9d6262",
          V: "#07f8f8",
          B: "#f30c0c",
          X: "#7c8383",
          Z: "#5ba4a4"
      };
      __webpack_require__(293);
      const pid = {
          init: function() {}
      }
        , gaps = ["", " ", "-", "_", "."];
      pid.run = function(letter, pos, conservation) {
          if (!conservation || 1 !== conservation.progress || pos > conservation.map.length || gaps.includes(letter))
              return "#ffffff";
          var cons = conservation.map[pos][letter] || 0;
          return .8 < cons ? "#6464ff" : .6 < cons ? "#9da5ff" : .4 < cons ? "#cccccc" : "#ffffff"
      }
      ,
      pid.map = {
          "> 0.8": "#6464ff",
          "> 0.6": "#9da5ff",
          "> 0.4": "#cccccc",
          "> 0": "#ffffff"
      };
      const staticSchemes = {
          aliphatic: {
              A: "#00639a",
              V: "#00639a",
              L: "#00639a",
              I: "#00639a",
              P: "#00639a"
          },
          aromatic: {
              F: "#00639a",
              W: "#00639a",
              Y: "#00639a",
              H: "#00639a"
          },
          buried: buried,
          buried_index: buried,
          charged: {
              D: "#00639a",
              E: "#00639a",
              K: "#00639a",
              R: "#00639a",
              H: "#00639a"
          },
          cinema: {
              A: "#BBBBBB",
              B: "grey",
              C: "yellow",
              D: "red",
              E: "red",
              F: "magenta",
              G: "brown",
              H: "#00FFFF",
              I: "#BBBBBB",
              J: "#fff",
              K: "#00FFFF",
              L: "#BBBBBB",
              M: "#BBBBBB",
              N: "green",
              O: "#fff",
              P: "brown",
              Q: "green",
              R: "#00FFFF",
              S: "green",
              T: "green",
              U: "#fff",
              V: "#BBBBBB",
              W: "magenta",
              X: "grey",
              Y: "magenta",
              Z: "grey",
              Gap: "grey"
          },
          clustal2: {
              A: "#80a0f0",
              R: "#f01505",
              N: "#00ff00",
              D: "#c048c0",
              C: "#f08080",
              Q: "#00ff00",
              E: "#c048c0",
              G: "#f09048",
              H: "#15a4a4",
              I: "#80a0f0",
              L: "#80a0f0",
              K: "#f01505",
              M: "#80a0f0",
              F: "#80a0f0",
              P: "#ffff00",
              S: "#00ff00",
              T: "#00ff00",
              W: "#80a0f0",
              Y: "#15a4a4",
              V: "#80a0f0",
              B: "#fff",
              X: "#fff",
              Z: "#fff"
          },
          clustal: {
              A: "orange",
              B: "#fff",
              C: "green",
              D: "red",
              E: "red",
              F: "blue",
              G: "orange",
              H: "red",
              I: "green",
              J: "#fff",
              K: "red",
              L: "green",
              M: "green",
              N: "#fff",
              O: "#fff",
              P: "orange",
              Q: "#fff",
              R: "red",
              S: "orange",
              T: "orange",
              U: "#fff",
              V: "green",
              W: "blue",
              X: "#fff",
              Y: "blue",
              Z: "#fff",
              Gap: "#fff"
          },
          helix: helix,
          helix_propensity: helix,
          hydro: {
              A: "#ad0052",
              B: "#0c00f3",
              C: "#c2003d",
              D: "#0c00f3",
              E: "#0c00f3",
              F: "#cb0034",
              G: "#6a0095",
              H: "#1500ea",
              I: "#ff0000",
              J: "#fff",
              K: "#0000ff",
              L: "#ea0015",
              M: "#b0004f",
              N: "#0c00f3",
              O: "#fff",
              P: "#4600b9",
              Q: "#0c00f3",
              R: "#0000ff",
              S: "#5e00a1",
              T: "#61009e",
              U: "#fff",
              V: "#f60009",
              W: "#5b00a4",
              X: "#680097",
              Y: "#4f00b0",
              Z: "#0c00f3"
          },
          lesk: {
              A: " orange",
              B: " #fff",
              C: " green",
              D: " red",
              E: " red",
              F: " green",
              G: " orange",
              H: " magenta",
              I: " green",
              J: " #fff",
              K: " red",
              L: " green",
              M: " green",
              N: " magenta",
              O: " #fff",
              P: " green",
              Q: " magenta",
              R: " red",
              S: " orange",
              T: " orange",
              U: " #fff",
              V: " green",
              W: " green",
              X: " #fff",
              Y: " green",
              Z: " #fff",
              Gap: " #fff"
          },
          mae: {
              A: " #77dd88",
              B: " #fff",
              C: " #99ee66",
              D: " #55bb33",
              E: " #55bb33",
              F: " #9999ff",
              G: " #77dd88",
              H: " #5555ff",
              I: " #66bbff",
              J: " #fff",
              K: " #ffcc77",
              L: " #66bbff",
              M: " #66bbff",
              N: " #55bb33",
              O: " #fff",
              P: " #eeaaaa",
              Q: " #55bb33",
              R: " #ffcc77",
              S: " #ff4455",
              T: " #ff4455",
              U: " #fff",
              V: " #66bbff",
              W: " #9999ff",
              X: " #fff",
              Y: " #9999ff",
              Z: " #fff",
              Gap: " #fff"
          },
          negative: {
              D: "#00639a",
              E: "#00639a"
          },
          nucleotide: {
              A: " #64F73F",
              C: " #FFB340",
              G: " #EB413C",
              T: " #3C88EE",
              U: " #3C88EE"
          },
          polar: {
              R: "#00639a",
              N: "#00639a",
              D: "#00639a",
              Q: "#00639a",
              E: "#00639a",
              G: "#00639a",
              H: "#00639a",
              K: "#00639a",
              S: "#00639a",
              T: "#00639a",
              Y: "#00639a"
          },
          positive: {
              K: "#00639a",
              R: "#00639a",
              H: "#00639a"
          },
          purine: purine,
          purine_pyrimidine: purine,
          serine_threonine: {
              S: "#00639a",
              T: "#00639a"
          },
          strand: strand,
          strand_propensity: strand,
          taylor: {
              A: "#ccff00",
              R: "#0000ff",
              N: "#cc00ff",
              D: "#ff0000",
              C: "#ffff00",
              Q: "#ff00cc",
              E: "#ff0066",
              G: "#ff9900",
              H: "#0066ff",
              I: "#66ff00",
              L: "#33ff00",
              K: "#6600ff",
              M: "#00ff00",
              F: "#00ff66",
              P: "#ffcc00",
              S: "#ff3300",
              T: "#ff6600",
              W: "#00ccff",
              Y: "#00ffcc",
              V: "#99ff00",
              B: "#fff",
              X: "#fff",
              Z: "#fff"
          },
          turn: turn,
          turn_propensity: turn,
          zappo: {
              A: "#ffafaf",
              R: "#6464ff",
              N: "#00ff00",
              D: "#ff0000",
              C: "#ffff00",
              Q: "#00ff00",
              E: "#ff0000",
              G: "#ff00ff",
              H: "#6464ff",
              I: "#ffafaf",
              L: "#ffafaf",
              K: "#6464ff",
              M: "#ffafaf",
              F: "#ffc800",
              P: "#ff00ff",
              S: "#00ff00",
              T: "#00ff00",
              W: "#ffc800",
              Y: "#ffc800",
              V: "#ffafaf",
              B: "#fff",
              X: "#fff",
              Z: "#fff"
          }
      }
        , dynSchemes = {
          conservation: pid
      }
        , Colors = function(opt) {
          this.maps = clone(staticSchemes),
          this.dyn = clone(dynSchemes),
          this.opt = opt
      };
      function clone(obj) {
          if (null == obj || "object" != typeof obj)
              return obj;
          var copy = obj.constructor();
          for (var attr in obj)
              obj.hasOwnProperty(attr) && (copy[attr] = obj[attr]);
          return copy
      }
      Colors.getScheme = function(scheme) {
          return staticSchemes[scheme]
      }
      ,
      Colors.prototype.getScheme = function(scheme) {
          var color = this.maps[scheme];
          return void 0 === color && (color = {},
          void 0 !== this.dyn[scheme]) ? new DynSchemeClass(this.dyn[scheme],this.opt) : new StaticSchemeClass(color)
      }
      ,
      Colors.prototype.addStaticScheme = function(name, scheme) {
          this.maps[name] = scheme
      }
      ,
      Colors.prototype.addDynScheme = function(name, scheme) {
          this.dyn[name] = scheme
      }
      ;
      __webpack_exports__.a = Colors
  },
  229: function(module, exports, __webpack_require__) {
      module.exports = function() {
          return new Worker(__webpack_require__.p + "c23331654b96503196f8.worker.js")
      }
  },
  230: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__(30);
      var prop_types = __webpack_require__(1)
        , prop_types_default = __webpack_require__.n(prop_types)
        , redux = __webpack_require__(138)
        , merge = __webpack_require__(925)
        , PropTypes = __webpack_require__(57)
        , actions = (__webpack_require__(15),
      __webpack_require__(21),
      __webpack_require__(5),
      __webpack_require__(8));
      var reduce = __webpack_require__(927);
      var reducers_calculateSequencesState = (prevState,sequences)=>({
          raw: sequences,
          length: sequences.length,
          maxLength: Object(reduce.a)(sequences, (m,e)=>Math.max(m, e.sequence.length), 0)
      })
        , ColorScheme = __webpack_require__(112);
      function ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function _objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                  _defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }
      function _defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      function checkColorScheme(state) {
          var _state$colorScheme;
          Object(ColorScheme.b)(state.colorScheme) || (state.colorScheme = new ColorScheme.a(state.colorScheme)),
          "function" == typeof (null == state || null === (_state$colorScheme = state.colorScheme) || void 0 === _state$colorScheme ? void 0 : _state$colorScheme.updateConservation) && state.colorScheme.updateConservation(state.conservation)
      }
      var reducers = ((reducer,statReducers)=>{
          const keys = Object.keys(statReducers);
          return function(prevState={}, action) {
              const state = reducer(prevState, action);
              if (prevState !== state)
                  for (let i = 0; i < keys.length; i++) {
                      const key = keys[i]
                        , nextStateForKey = statReducers[key](state[key], action, state);
                      if (void 0 === nextStateForKey)
                          throw new Error("A reducer can't return 'undefined'");
                      state[key] = nextStateForKey
                  }
              return state
          }
      }
      )(function combineReducers(reducers) {
          const keys = Object.keys(reducers);
          return function(state={}, action) {
              const nextState = {};
              let hasChanged = !1;
              for (let i = 0; i < keys.length; i++) {
                  const key = keys[i]
                    , nextStateForKey = reducers[key](state[key], action);
                  if (void 0 === nextStateForKey)
                      throw new Error("A reducer can't return 'undefined'");
                  nextState[key] = nextStateForKey,
                  hasChanged = hasChanged || nextStateForKey !== state[key]
              }
              return hasChanged ? nextState : state
          }
      }({
          props: (state={},{type: type, payload: payload})=>{
              switch (type) {
              case actions.b.updateProps.key:
                  return state = _objectSpread({}, state, {}, payload),
                  "colorScheme"in payload && checkColorScheme(state),
                  state;
              case actions.b.updateProp.key:
                  const key = payload.key;
                  return state = _objectSpread({}, state, {
                      [key]: payload.value
                  }),
                  "colorScheme" === key && checkColorScheme(state),
                  state;
              case actions.b.updateConservation.key:
                  const progress = payload.progress
                    , conservation = payload.conservation;
                  return state.conservation = {
                      progress: progress,
                      map: conservation
                  },
                  state.colorScheme.updateConservation(state.conservation),
                  state;
              default:
                  return void 0 !== state.colorScheme && checkColorScheme(state),
                  state
              }
          }
          ,
          sequences: function handleActions(handlers, initialState) {
              return function(state=initialState, {type: type, payload: payload}) {
                  return handlers.hasOwnProperty(type) ? handlers[type](state, payload) : state
              }
          }({
              [actions.b.updateSequences]: reducers_calculateSequencesState
          }, [])
      }), {
          sequenceStats: (prevState={
              currentViewSequence: 0,
              currentViewSequencePosition: 0
          },action,state)=>{
              switch (action.type) {
              case actions.b.updateProp.key:
              case actions.b.updateProps.key:
              case actions.b.updateSequences.key:
                  if (state.props && state.props.tileHeight && state.props.tileWidth && state.sequences) {
                      return {
                          nrXTiles: Math.ceil(state.props.width / state.props.tileWidth) + 1,
                          nrYTiles: Math.ceil(state.props.height / state.props.tileHeight) + 1,
                          fullWidth: state.props.tileWidth * state.sequences.maxLength,
                          fullHeight: state.props.tileHeight * state.sequences.length
                      }
                  }
              }
              return prevState
          }
      })
        , debug = __webpack_require__(158);
      function _objectWithoutProperties(source, excluded) {
          if (null == source)
              return {};
          var key, i, target = function _objectWithoutPropertiesLoose(source, excluded) {
              if (null == source)
                  return {};
              var key, i, target = {}, sourceKeys = Object.keys(source);
              for (i = 0; i < sourceKeys.length; i++)
                  key = sourceKeys[i],
                  0 <= excluded.indexOf(key) || (target[key] = source[key]);
              return target
          }(source, excluded);
          if (Object.getOwnPropertySymbols) {
              var sourceSymbolKeys = Object.getOwnPropertySymbols(source);
              for (i = 0; i < sourceSymbolKeys.length; i++)
                  key = sourceSymbolKeys[i],
                  0 <= excluded.indexOf(key) || Object.prototype.propertyIsEnumerable.call(source, key) && (target[key] = source[key])
          }
          return target
      }
      __webpack_exports__.a = props=>{
          prop_types_default.a.checkPropTypes(PropTypes.b, props, "prop", "MSAViewer");
          const propsWithDefaultValues = Object(merge.a)({}, PropTypes.e, props)
            , sequences = propsWithDefaultValues.sequences
            , otherProps = (propsWithDefaultValues.position,
          _objectWithoutProperties(propsWithDefaultValues, ["sequences", "position"]))
            , store = Object(redux.b)(reducers, debug.a && window.__REDUX_DEVTOOLS_EXTENSION__ && window.__REDUX_DEVTOOLS_EXTENSION__());
          return store.dispatch(Object(actions.c)(otherProps)),
          store.dispatch(Object(actions.d)(sequences)),
          store
      }
  },
  246: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__(33);
      var react__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(0)
        , react__WEBPACK_IMPORTED_MODULE_1___default = __webpack_require__.n(react__WEBPACK_IMPORTED_MODULE_1__)
        , prop_types__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(1)
        , prop_types__WEBPACK_IMPORTED_MODULE_2___default = __webpack_require__.n(prop_types__WEBPACK_IMPORTED_MODULE_2__);
      class ListComponent extends react__WEBPACK_IMPORTED_MODULE_1__.PureComponent {
          renderTile(i) {
              const TileComponent = this.props.tileComponent;
              if (i in this.props.componentCache)
                  return this.props.componentCache[i];
              {
                  const el = react__WEBPACK_IMPORTED_MODULE_1___default.a.createElement(TileComponent, {
                      key: i,
                      index: i
                  });
                  return this.props.componentCache[i] = el,
                  el
              }
          }
          render() {
              const elements = [];
              for (let i = this.props.startTile; i < this.props.endTile; i++)
                  elements.push(this.renderTile(i));
              return 0 === elements.length && console.warn("The TileComponent rendered returned 0 elements from ".concat(this.props.startTile, " to ").concat(this.props.endTile)),
              react__WEBPACK_IMPORTED_MODULE_1___default.a.createElement("div", null, elements)
          }
      }
      ListComponent.displayName = "ListComponent",
      ListComponent.propTypes = {
          startTile: prop_types__WEBPACK_IMPORTED_MODULE_2___default.a.number.isRequired,
          endTile: prop_types__WEBPACK_IMPORTED_MODULE_2___default.a.number.isRequired,
          componentCache: prop_types__WEBPACK_IMPORTED_MODULE_2___default.a.func.isRequired
      },
      ListComponent.__docgenInfo = {
          description: "Renders a list of tiles, but caches already seen components.",
          methods: [{
              name: "renderTile",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "i",
                  type: null
              }],
              returns: null
          }],
          displayName: "ListComponent",
          props: {
              startTile: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              endTile: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              componentCache: {
                  type: {
                      name: "func"
                  },
                  required: !0,
                  description: ""
              }
          }
      },
      __webpack_exports__.a = ListComponent,
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/ListComponent.js"] = {
          name: "ListComponent",
          docgenInfo: ListComponent.__docgenInfo,
          path: "src/components/ListComponent.js"
      })
  },
  249: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__(15),
      __webpack_require__(30),
      __webpack_require__(7),
      __webpack_require__(21),
      __webpack_require__(5);
      var react__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(0)
        , react__WEBPACK_IMPORTED_MODULE_5___default = __webpack_require__.n(react__WEBPACK_IMPORTED_MODULE_5__)
        , prop_types__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(1)
        , prop_types__WEBPACK_IMPORTED_MODULE_6___default = __webpack_require__.n(prop_types__WEBPACK_IMPORTED_MODULE_6__)
        , _store_withPositionStore__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(37)
        , _ListComponent__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(246);
      function _extends() {
          return (_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      function ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function _defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      function _objectWithoutProperties(source, excluded) {
          if (null == source)
              return {};
          var key, i, target = function _objectWithoutPropertiesLoose(source, excluded) {
              if (null == source)
                  return {};
              var key, i, target = {}, sourceKeys = Object.keys(source);
              for (i = 0; i < sourceKeys.length; i++)
                  key = sourceKeys[i],
                  0 <= excluded.indexOf(key) || (target[key] = source[key]);
              return target
          }(source, excluded);
          if (Object.getOwnPropertySymbols) {
              var sourceSymbolKeys = Object.getOwnPropertySymbols(source);
              for (i = 0; i < sourceSymbolKeys.length; i++)
                  key = sourceSymbolKeys[i],
                  0 <= excluded.indexOf(key) || Object.prototype.propertyIsEnumerable.call(source, key) && (target[key] = source[key])
          }
          return target
      }
      class XBarComponent extends react__WEBPACK_IMPORTED_MODULE_5__.PureComponent {
          constructor(props) {
              super(props),
              this.el = react__WEBPACK_IMPORTED_MODULE_5___default.a.createRef()
          }
          render() {
              const _this$props = this.props
                , tileWidth = _this$props.tileWidth
                , sequences = _this$props.sequences
                , width = _this$props.width
                , tileComponent = (_this$props.cacheElements,
              _this$props.tileComponent)
                , otherProps = (_this$props.nrXTiles,
              _this$props.maxLength,
              _this$props.position,
              _this$props.positionDispatch,
              _this$props.componentCache,
              _objectWithoutProperties(_this$props, ["tileWidth", "sequences", "width", "cacheElements", "tileComponent", "nrXTiles", "maxLength", "position", "positionDispatch", "componentCache"]))
                , containerStyle = function _objectSpread(target) {
                  for (var source, i = 1; i < arguments.length; i++)
                      source = null != arguments[i] ? arguments[i] : {},
                      i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                          _defineProperty(target, key, source[key])
                      }
                      )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                          Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
                      }
                      ));
                  return target
              }({}, this.props.style, {
                  height: this.props.height
              })
                , startTile = Math.max(0, this.props.position.currentViewSequencePosition - this.props.cacheElements)
                , endTile = Math.min(this.props.maxLength, startTile + this.props.nrXTiles + 2 * this.props.cacheElements)
                , maxWidth = this.props.width + 2 * this.props.cacheElements * this.props.tileWidth;
              return this.props.position.lastStartXTile = startTile,
              this.props.position.lastCurrentViewSequencePosition = this.props.position.currentViewSequencePosition,
              react__WEBPACK_IMPORTED_MODULE_5___default.a.createElement("div", _extends({
                  style: containerStyle
              }, otherProps), react__WEBPACK_IMPORTED_MODULE_5___default.a.createElement("div", {
                  style: {
                      width: width,
                      overflow: "hidden",
                      position: "relative",
                      whiteSpace: "nowrap"
                  },
                  ref: this.el
              }, react__WEBPACK_IMPORTED_MODULE_5___default.a.createElement(_ListComponent__WEBPACK_IMPORTED_MODULE_8__.a, _extends({}, {
                  tileWidth: tileWidth,
                  sequences: sequences,
                  tileComponent: tileComponent
              }, {
                  componentCache: this.props.componentCache,
                  startTile: startTile,
                  endTile: endTile,
                  maxWidth: maxWidth
              }))))
          }
      }
      XBarComponent.displayName = "XBarComponent",
      XBarComponent.propTypes = {
          tileComponent: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.oneOfType([prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.func, prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.object]).isRequired,
          cacheElements: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.number.isRequired,
          tileWidth: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.number.isRequired,
          nrXTiles: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.number.isRequired,
          maxLength: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.number.isRequired,
          componentCache: prop_types__WEBPACK_IMPORTED_MODULE_6___default.a.func.isRequired
      },
      XBarComponent.__docgenInfo = {
          description: "Displays the sequence names with an arbitrary Marker component",
          methods: [],
          displayName: "XBarComponent",
          props: {
              tileComponent: {
                  type: {
                      name: "union",
                      value: [{
                          name: "func"
                      }, {
                          name: "object"
                      }]
                  },
                  required: !0,
                  description: "Tile to render."
              },
              cacheElements: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              tileWidth: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              nrXTiles: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              maxLength: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              componentCache: {
                  type: {
                      name: "func"
                  },
                  required: !0,
                  description: ""
              }
          }
      },
      __webpack_exports__.a = Object(_store_withPositionStore__WEBPACK_IMPORTED_MODULE_7__.a)(XBarComponent, {
          withX: !0
      }),
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/xBars/xBar.js"] = {
          name: "XBarComponent",
          docgenInfo: XBarComponent.__docgenInfo,
          path: "src/components/xBars/xBar.js"
      })
  },
  254: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__(15),
      __webpack_require__(30),
      __webpack_require__(7),
      __webpack_require__(21),
      __webpack_require__(5);
      var react = __webpack_require__(0)
        , react_default = __webpack_require__.n(react)
        , prop_types = __webpack_require__(1)
        , prop_types_default = __webpack_require__.n(prop_types)
        , partialRight = __webpack_require__(915)
        , pick = __webpack_require__(319)
        , connect = __webpack_require__(16)
        , shallowSelect = __webpack_require__(76)
        , autobind = __webpack_require__(77)
        , withPositionStore = __webpack_require__(37)
        , ListComponent = __webpack_require__(246);
      function _objectWithoutProperties(source, excluded) {
          if (null == source)
              return {};
          var key, i, target = function _objectWithoutPropertiesLoose(source, excluded) {
              if (null == source)
                  return {};
              var key, i, target = {}, sourceKeys = Object.keys(source);
              for (i = 0; i < sourceKeys.length; i++)
                  key = sourceKeys[i],
                  0 <= excluded.indexOf(key) || (target[key] = source[key]);
              return target
          }(source, excluded);
          if (Object.getOwnPropertySymbols) {
              var sourceSymbolKeys = Object.getOwnPropertySymbols(source);
              for (i = 0; i < sourceSymbolKeys.length; i++)
                  key = sourceSymbolKeys[i],
                  0 <= excluded.indexOf(key) || Object.prototype.propertyIsEnumerable.call(source, key) && (target[key] = source[key])
          }
          return target
      }
      class yBar_YBarComponent extends react.PureComponent {
          constructor(props) {
              super(props),
              this.el = react_default.a.createRef()
          }
          render() {
              const _this$props = this.props
                , height = (_this$props.yPosOffset,
              _this$props.tileHeight,
              _this$props.currentViewSequence,
              _this$props.sequences,
              _this$props.height)
                , otherProps = (_this$props.cacheElements,
              _this$props.tileComponent,
              _this$props.nrYTiles,
              _this$props.position,
              _this$props.positionDispatch,
              _this$props.componentCache,
              _objectWithoutProperties(_this$props, ["yPosOffset", "tileHeight", "currentViewSequence", "sequences", "height", "cacheElements", "tileComponent", "nrYTiles", "position", "positionDispatch", "componentCache"]))
                , startTile = Math.max(0, this.props.position.currentViewSequence - this.props.cacheElements)
                , endTile = Math.min(this.props.sequences.length, startTile + Math.ceil(height / this.props.tileHeight) + 2 * this.props.cacheElements);
              return this.props.position.lastCurrentViewSequence = this.props.position.currentViewSequence,
              this.props.position.lastStartYTile = startTile,
              react_default.a.createElement("div", otherProps, react_default.a.createElement("div", {
                  style: {
                      height: height,
                      overflow: "hidden",
                      position: "relative",
                      whiteSpace: "nowrap"
                  },
                  ref: this.el
              }, react_default.a.createElement(ListComponent.a, {
                  componentCache: this.props.componentCache,
                  startTile: startTile,
                  endTile: endTile,
                  tileComponent: this.props.tileComponent
              })))
          }
      }
      yBar_YBarComponent.displayName = "YBarComponent",
      yBar_YBarComponent.propTypes = {
          tileComponent: prop_types_default.a.oneOfType([prop_types_default.a.func, prop_types_default.a.object]).isRequired,
          cacheElements: prop_types_default.a.number.isRequired,
          tileHeight: prop_types_default.a.number.isRequired,
          nrYTiles: prop_types_default.a.number.isRequired,
          componentCache: prop_types_default.a.func.isRequired
      },
      yBar_YBarComponent.__docgenInfo = {
          description: "Displays the sequence names with an arbitrary Marker component",
          methods: [],
          displayName: "YBarComponent",
          props: {
              tileComponent: {
                  type: {
                      name: "union",
                      value: [{
                          name: "func"
                      }, {
                          name: "object"
                      }]
                  },
                  required: !0,
                  description: "Tile to render."
              },
              cacheElements: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              tileHeight: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              nrYTiles: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              componentCache: {
                  type: {
                      name: "func"
                  },
                  required: !0,
                  description: ""
              }
          }
      };
      var yBar = Object(withPositionStore.a)(yBar_YBarComponent, {
          withY: !0
      });
      function _extends() {
          return (_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      function ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function _objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                  _defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }
      function _defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      function Labels_objectWithoutProperties(source, excluded) {
          if (null == source)
              return {};
          var key, i, target = function Labels_objectWithoutPropertiesLoose(source, excluded) {
              if (null == source)
                  return {};
              var key, i, target = {}, sourceKeys = Object.keys(source);
              for (i = 0; i < sourceKeys.length; i++)
                  key = sourceKeys[i],
                  0 <= excluded.indexOf(key) || (target[key] = source[key]);
              return target
          }(source, excluded);
          if (Object.getOwnPropertySymbols) {
              var sourceSymbolKeys = Object.getOwnPropertySymbols(source);
              for (i = 0; i < sourceSymbolKeys.length; i++)
                  key = sourceSymbolKeys[i],
                  0 <= excluded.indexOf(key) || Object.prototype.propertyIsEnumerable.call(source, key) && (target[key] = source[key])
          }
          return target
      }
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/yBars/yBar.js"] = {
          name: "YBarComponent",
          docgenInfo: yBar_YBarComponent.__docgenInfo,
          path: "src/components/yBars/yBar.js"
      });
      class Labels_HTMLLabelsComponent extends react.PureComponent {
          constructor(props) {
              super(props),
              Object(autobind.a)(this, "createLabel"),
              this.label = Object(shallowSelect.a)(Object(partialRight.a)(pick.a, this.constructor.labelProps), this.createLabel)
          }
          createLabel(props) {
              return this.cache = function() {}
              ,
              function createLabel({sequences: sequences, tileHeight: tileHeight, labelComponent: labelComponent, labelStyle: labelStyle, labelAttributes: labelAttributes}) {
                  class Label extends react.PureComponent {
                      render() {
                          const _this$props = this.props
                            , index = _this$props.index
                            , otherProps = Labels_objectWithoutProperties(_this$props, ["index"]);
                          if (labelComponent)
                              return react_default.a.createElement(labelComponent, {
                                  sequence: sequences[index],
                                  index: index
                              });
                          {
                              otherProps.style = _objectSpread({}, this.props.style, {
                                  height: tileHeight
                              }, labelStyle);
                              const attributes = _objectSpread({}, otherProps, {}, labelAttributes);
                              return react_default.a.createElement("div", attributes, sequences[index].name)
                          }
                      }
                  }
                  return Label.displayName = "Label",
                  Label
              }(props)
          }
          render() {
              const _this$props2 = this.props
                , cacheElements = _this$props2.cacheElements
                , otherProps = (_this$props2.dispatch,
              _this$props2.labelComponent,
              _this$props2.labelStyle,
              _this$props2.labelAttributes,
              Labels_objectWithoutProperties(_this$props2, ["cacheElements", "dispatch", "labelComponent", "labelStyle", "labelAttributes"]));
              return react_default.a.createElement(yBar, _extends({
                  tileComponent: this.label(this.props),
                  cacheElements: cacheElements,
                  componentCache: this.cache
              }, otherProps))
          }
      }
      Labels_HTMLLabelsComponent.displayName = "HTMLLabelsComponent",
      _defineProperty(Labels_HTMLLabelsComponent, "labelProps", ["sequences", "tileHeight", "labelComponent", "labelStyle", "labelAttributes"]),
      Labels_HTMLLabelsComponent.defaultProps = {
          cacheElements: 10,
          labelStyle: {}
      },
      Labels_HTMLLabelsComponent.propTypes = {
          font: prop_types_default.a.string,
          labelComponent: prop_types_default.a.oneOfType([prop_types_default.a.object, prop_types_default.a.func]),
          style: prop_types_default.a.object,
          labelStyle: prop_types_default.a.object,
          labelAttributes: prop_types_default.a.object
      };
      Labels_HTMLLabelsComponent.__docgenInfo = {
          description: "Displays the sequence names.",
          methods: [{
              name: "createLabel",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "props",
                  type: null
              }],
              returns: null
          }],
          displayName: "HTMLLabelsComponent",
          props: {
              cacheElements: {
                  defaultValue: {
                      value: "10",
                      computed: !1
                  },
                  required: !1
              },
              labelStyle: {
                  defaultValue: {
                      value: "{}",
                      computed: !1
                  },
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Inline styles to apply to each label."
              },
              font: {
                  type: {
                      name: "string"
                  },
                  required: !1,
                  description: "Font of the sequence labels, e.g. `20px Arial`"
              },
              labelComponent: {
                  type: {
                      name: "union",
                      value: [{
                          name: "object"
                      }, {
                          name: "func"
                      }]
                  },
                  required: !1,
                  description: "Component to create labels from."
              },
              style: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Inline styles to apply to the Label component"
              },
              labelAttributes: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Attributes to apply to each label."
              }
          }
      };
      __webpack_exports__.a = Object(connect.a)(state=>({
          height: state.props.height,
          tileHeight: state.props.tileHeight,
          sequences: state.sequences.raw,
          nrYTiles: state.sequenceStats.nrYTiles
      }))(Labels_HTMLLabelsComponent);
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/yBars/Labels.js"] = {
          name: "HTMLLabelsComponent",
          docgenInfo: Labels_HTMLLabelsComponent.__docgenInfo,
          path: "src/components/yBars/Labels.js"
      })
  },
  255: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__(15),
      __webpack_require__(30),
      __webpack_require__(35),
      __webpack_require__(7),
      __webpack_require__(21),
      __webpack_require__(5);
      var react = __webpack_require__(0)
        , react_default = __webpack_require__.n(react)
        , prop_types = __webpack_require__(1)
        , prop_types_default = __webpack_require__.n(prop_types)
        , pick = __webpack_require__(319)
        , es = __webpack_require__(185)
        , xBar = __webpack_require__(249)
        , connect = __webpack_require__(16)
        , shallowSelect = __webpack_require__(76)
        , autobind = __webpack_require__(77)
        , assignIn = (__webpack_require__(459),
      __webpack_require__(919))
        , forEach = __webpack_require__(920)
        , mapValues = __webpack_require__(921)
        , map = __webpack_require__(931)
        , reduce = __webpack_require__(927)
        , reject = __webpack_require__(926)
        , max = __webpack_require__(929);
      const stat = function(seqs, opts) {
          if (!this || this.constructor !== stat)
              return new stat(seqs);
          if (void 0 === seqs || "string" == typeof seqs)
              throw new TypeError("you need to give the seq stat an array");
          this.resetSeqs(seqs),
          this.alphabetSize = 4,
          this._useBackground = !1,
          this.useGaps = !1,
          this.ignoredChars = ["-", "*"],
          Object(assignIn.a)(this, opts)
      };
      stat.prototype.addSeq = function(seq) {
          this.seqs.push(seq),
          this._reset()
      }
      ,
      stat.prototype.removeSeq = function(seq) {
          "number" == typeof seq ? this.seqs.splice(seq, 1) : Object(forEach.a)(this.seqs, function(s, i) {
              seq === s && this.seqs.splice(i, 1)
          }
          .bind(this)),
          this._reset()
      }
      ,
      stat.prototype.addSeqs = function(seqs) {
          seqs.forEach(function(seq) {
              this.addSeq(seq)
          }
          .bind(this))
      }
      ,
      stat.prototype.resetSeqs = function(seqs) {
          if (this.seqs = [],
          !seqs instanceof Array || "at"in seqs) {
              this.mseqs = seqs;
              var mSeqsPluck = function() {
                  var seqArr = this.mseqs.pluck("seq");
                  this.resetSeqs(seqArr)
              };
              seqs.on("add change reset ", mSeqsPluck, this),
              mSeqsPluck.call(this)
          } else
              this.addSeqs(seqs),
              this._reset()
      }
      ;
      var calcValues = ["consensus", "frequency", "maxLength", "ic", "gaps"];
      stat.prototype._reset = function() {
          for (var i = 0; i < calcValues.length; i++)
              this["_" + calcValues[i]] = void 0;
          this._identity = void 0,
          this._background = void 0
      }
      ,
      stat.prototype.setBackground = function(b) {
          this._useBackground = b,
          this._reset()
      }
      ,
      stat.prototype.useBackground = function() {
          this.setBackground(!0)
      }
      ,
      stat.prototype.setDNA = function() {
          this.alphabetSize = 4
      }
      ,
      stat.prototype.setProtein = function() {
          this.alphabetSize = 20
      }
      ,
      calcValues.forEach((function(key) {
          stat.prototype[key] = function() {
              return void 0 === this["_" + key] && (this["_" + key] = this[key + "Calc"]()),
              this["_" + key]
          }
      }
      )),
      stat.prototype.identity = function(seq) {
          var ident;
          return (void 0 === this._identity || seq) && (ident = this.identityCalc(seq),
          this._identity = void 0),
          this._identity || ident
      }
      ,
      stat.prototype.background = function() {
          return void 0 !== this.bg ? this.bg : (void 0 === this._background && this.backgroundCalc(),
          this._background)
      }
      ,
      stat.prototype.frequencyCalc = function(opts) {
          var occs, totalPerPos;
          occs = Array(this.maxLength()),
          totalPerPos = Array(this.seqs.length);
          var ignoredChars = this.ignoredChars;
          return void 0 !== opts && opts.all && (ignoredChars = []),
          Object(forEach.a)(this.seqs, (function(el) {
              Object(forEach.a)(el, (function(c, pos) {
                  0 <= ignoredChars.indexOf(c) || (void 0 === occs[pos] && (occs[pos] = {}),
                  void 0 === occs[pos][c] && (occs[pos][c] = 0),
                  occs[pos][c]++,
                  void 0 === totalPerPos[pos] && (totalPerPos[pos] = 0),
                  totalPerPos[pos]++)
              }
              ))
          }
          )),
          Object(forEach.a)(occs, (function(el, pos) {
              return Object(forEach.a)(el, (function(val, c) {
                  return occs[pos][c] = val / totalPerPos[pos]
              }
              ))
          }
          )),
          this._frequency = occs,
          occs
      }
      ,
      stat.prototype.backgroundCalc = function() {
          var occ = {}
            , total = 0;
          return Object(forEach.a)(this.seqs, (function(el) {
              Object(forEach.a)(el, (function(c) {
                  return void 0 === occ[c] && (occ[c] = 0),
                  occ[c]++,
                  total++
              }
              ))
          }
          )),
          occ = Object(mapValues.a)(occ, (function(val) {
              return val / total
          }
          )),
          this._background = occ,
          occ
      }
      ,
      stat.prototype.icCalc = function() {
          var f = this.frequency();
          if (this._useBackground)
              var b = this.background();
          var ignoredChars = this.ignoredChars
            , useBackground = this._useBackground
            , ic = Object(map.a)(f, (function(el) {
              return Object(reduce.a)(el, (function(memo, val, c) {
                  return 0 <= ignoredChars.indexOf(c) ? memo : (useBackground && (val /= b[c]),
                  memo - val * (Math.log(val) / Math.log(2)))
              }
              ), 0)
          }
          ));
          return this._ic = ic,
          ic
      }
      ,
      stat.prototype.conservation = function(alphabetSize) {
          var ic = this.ic()
            , gaps = this.gaps()
            , self = this;
          alphabetSize = alphabetSize || this.alphabetSize;
          var icMax = Math.log(alphabetSize) / Math.log(2)
            , i = 0;
          return Object(map.a)(ic, (function(el) {
              var ret = icMax - el;
              return self.useGaps && (ret *= 1 - gaps[i++]),
              ret
          }
          ))
      }
      ,
      stat.prototype.conservResidue = function(input) {
          var ic, alphabetSize = input ? input.alphabetSize : void 0, ignoredChars = this.ignoredChars;
          ic = void 0 !== input && input.scaled ? this.scale(this.conservation(alphabetSize)) : this.conservation(alphabetSize);
          var keys, f = this.frequency();
          return Object(map.a)(f, (function(el, i) {
              keys = Object(reject.a)(keys(el), (function(c) {
                  return 0 <= ignoredChars.indexOf(c)
              }
              ));
              var obj = {};
              return Object(forEach.a)(keys, (function(key) {
                  obj[key] = el[key] * ic[i]
              }
              )),
              obj
          }
          ))
      }
      ,
      stat.prototype.conservResidue2 = function(alphabetSize) {
          var f = this.frequency()
            , ic = this.conservation(alphabetSize)
            , b = this.background();
          return Object(map.a)(f, (function(el, i) {
              return Object(map.a)(el, (function(val) {
                  var sum = Object(reduce.a)(f[i], (function(memo, e) {
                      return memo + e / b[i]
                  }
                  ), 0);
                  return val / b[i] / sum * ic[i]
              }
              ), 0)
          }
          ))
      }
      ,
      stat.prototype.scale = function(ic, alphabetSize) {
          alphabetSize = alphabetSize || this.alphabetSize;
          var icMax = Math.log(alphabetSize) / Math.log(2);
          return Object(map.a)(ic, (function(el) {
              return el / icMax
          }
          ))
      }
      ,
      stat.prototype.maxLengthCalc = function() {
          return 0 === this.seqs.length ? 0 : Object(max.a)(this.seqs, (function(seq) {
              return seq.length
          }
          )).length
      }
      ,
      stat.prototype.consensusCalc = function() {
          var occs = Array(this.maxLength());
          return Object(forEach.a)(this.seqs, (function(el) {
              Object(forEach.a)(el, (function(c, pos) {
                  void 0 === occs[pos] && (occs[pos] = {}),
                  void 0 === occs[pos][c] && (occs[pos][c] = 0),
                  occs[pos][c]++
              }
              ))
          }
          )),
          this._consensus = Object(reduce.a)(occs, (function(memo, occ) {
              var keys = Object.keys(occ);
              return memo + Object(max.a)(keys, (function(key) {
                  return occ[key]
              }
              ))
          }
          ), ""),
          this._consensus
      }
      ,
      stat.prototype.identityCalc = function(compareSeq) {
          var consensus = compareSeq || this.consensus();
          return this._identity = this.seqs.map((function(seq) {
              for (var matches = 0, total = 0, i = 0; i < seq.length; i++)
                  "-" !== seq[i] && "-" !== consensus[i] && (total++,
                  seq[i] === consensus[i] && matches++);
              return matches / total
          }
          )),
          this._identity
      }
      ,
      stat.prototype.gapsCalc = function() {
          var mLength = this.maxLength();
          if (1 >= mLength || void 0 === mLength)
              return [];
          var occs = Array(this.maxLength());
          return Object(forEach.a)(this.seqs, (function(el) {
              Object(forEach.a)(el, (function(c, pos) {
                  void 0 === occs[pos] && (occs[pos] = {
                      g: 0,
                      t: 0
                  }),
                  c = "-" === c ? "g" : "t",
                  occs[pos][c]++
              }
              ))
          }
          )),
          this._gaps = Object(map.a)(occs, (function(el) {
              return el.g / (el.g + el.t)
          }
          )),
          this._gaps
      }
      ;
      var statSeqs = stat;
      function _extends() {
          return (_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      function ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function _defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      function _objectWithoutProperties(source, excluded) {
          if (null == source)
              return {};
          var key, i, target = function _objectWithoutPropertiesLoose(source, excluded) {
              if (null == source)
                  return {};
              var key, i, target = {}, sourceKeys = Object.keys(source);
              for (i = 0; i < sourceKeys.length; i++)
                  key = sourceKeys[i],
                  0 <= excluded.indexOf(key) || (target[key] = source[key]);
              return target
          }(source, excluded);
          if (Object.getOwnPropertySymbols) {
              var sourceSymbolKeys = Object.getOwnPropertySymbols(source);
              for (i = 0; i < sourceSymbolKeys.length; i++)
                  key = sourceSymbolKeys[i],
                  0 <= excluded.indexOf(key) || Object.prototype.propertyIsEnumerable.call(source, key) && (target[key] = source[key])
          }
          return target
      }
      class OverviewBar_HTMLOverviewBarComponent extends react.PureComponent {
          constructor(props) {
              super(props),
              this.cache = function() {}
              ,
              this.initializeColumnHeights(),
              Object(autobind.a)(this, "createBar"),
              this.bar = Object(shallowSelect.a)(s=>Object(pick.a)(s, this.constructor.barAttributes), this.columnHeights, this.createBar)
          }
          createBar(props, columnHeights) {
              return this.cache = function() {}
              ,
              function createBar({columnHeights: columnHeights, tileWidth: tileWidth, height: height, fillColor: fillColor, barStyle: barStyle, barAttributes: barAttributes}) {
                  class Bar extends react.PureComponent {
                      render() {
                          const _this$props = this.props
                            , index = _this$props.index
                            , otherProps = _objectWithoutProperties(_this$props, ["index"]);
                          return otherProps.style = {
                              height: Math.round(columnHeights[index] * height),
                              width: tileWidth,
                              display: "inline-block",
                              textAlign: "center",
                              backgroundColor: fillColor
                          },
                          react_default.a.createElement("div", otherProps)
                      }
                  }
                  return Bar.displayName = "Bar",
                  Bar
              }(function _objectSpread(target) {
                  for (var source, i = 1; i < arguments.length; i++)
                      source = null != arguments[i] ? arguments[i] : {},
                      i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                          _defineProperty(target, key, source[key])
                      }
                      )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                          Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
                      }
                      ));
                  return target
              }({}, props, {
                  columnHeights: columnHeights
              }))
          }
          initializeColumnHeights() {
              this.columnHeights = Object(es.a)(p=>p.sequences, p=>p.method, (sequences,method)=>{
                  const stats = statSeqs(sequences.map(e=>e.sequence));
                  let result;
                  switch (method) {
                  case "conservation":
                      result = stats.scale(stats.conservation());
                      break;
                  case "information-content":
                      result = stats.scale(stats.ic());
                      break;
                  default:
                      console.error(method + "is an invalid aggregation method for <OverviewBar />")
                  }
                  return result
              }
              ).bind(this)
          }
          render() {
              const _this$props2 = this.props
                , cacheElements = _this$props2.cacheElements
                , otherProps = (_this$props2.height,
              _this$props2.method,
              _this$props2.fillColor,
              _this$props2.dispatch,
              _this$props2.barStyle,
              _this$props2.barAttributes,
              _objectWithoutProperties(_this$props2, ["cacheElements", "height", "method", "fillColor", "dispatch", "barStyle", "barAttributes"]));
              return react_default.a.createElement(xBar.a, _extends({
                  tileComponent: this.bar(this.props),
                  cacheElements: cacheElements,
                  componentCache: this.cache
              }, otherProps))
          }
      }
      OverviewBar_HTMLOverviewBarComponent.displayName = "HTMLOverviewBarComponent",
      _defineProperty(OverviewBar_HTMLOverviewBarComponent, "barAttributes", ["tileWidth", "height", "fillColor", "barStyle", "barAttributes"]),
      OverviewBar_HTMLOverviewBarComponent.defaultProps = {
          height: 50,
          fillColor: "#999999",
          method: "conservation",
          cacheElements: 10
      },
      OverviewBar_HTMLOverviewBarComponent.propTypes = {
          method: prop_types_default.a.oneOf(["information-content", "conservation"]),
          height: prop_types_default.a.number,
          fillColor: prop_types_default.a.string,
          style: prop_types_default.a.object,
          barStyle: prop_types_default.a.object,
          barAttributes: prop_types_default.a.object
      };
      OverviewBar_HTMLOverviewBarComponent.__docgenInfo = {
          description: "Creates a small overview box of the sequences for a general overview.",
          methods: [{
              name: "createBar",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "props",
                  type: null
              }, {
                  name: "columnHeights",
                  type: null
              }],
              returns: null
          }, {
              name: "initializeColumnHeights",
              docblock: "Reduces the `props` object to column height by a `props.method`",
              modifiers: [],
              params: [],
              returns: null,
              description: "Reduces the `props` object to column height by a `props.method`"
          }],
          displayName: "HTMLOverviewBarComponent",
          props: {
              height: {
                  defaultValue: {
                      value: "50",
                      computed: !1
                  },
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Height of the OverviewBar (in pixels), e.g. `100`"
              },
              fillColor: {
                  defaultValue: {
                      value: '"#999999"',
                      computed: !1
                  },
                  type: {
                      name: "string"
                  },
                  required: !1,
                  description: "Fill color of the OverviewBar, e.g. `#999999`"
              },
              method: {
                  defaultValue: {
                      value: '"conservation"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"information-content"',
                          computed: !1
                      }, {
                          value: '"conservation"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: "Method to use for the OverviewBar:\n - `information-content`: Information entropy after Shannon of a column (scaled)\n - `conservation`: Conservation of a column (scaled)"
              },
              cacheElements: {
                  defaultValue: {
                      value: "10",
                      computed: !1
                  },
                  required: !1
              },
              style: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Inline styles to apply to the OverviewBar component"
              },
              barStyle: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Inline styles to apply to each bar."
              },
              barAttributes: {
                  type: {
                      name: "object"
                  },
                  required: !1,
                  description: "Attributes to apply to each bar."
              }
          }
      };
      __webpack_exports__.a = Object(connect.a)(state=>({
          sequences: state.sequences.raw,
          maxLength: state.sequences.maxLength,
          width: state.props.width,
          tileWidth: state.props.tileWidth,
          nrXTiles: state.sequenceStats.nrXTiles
      }))(OverviewBar_HTMLOverviewBarComponent);
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/xBars/OverviewBar.js"] = {
          name: "HTMLOverviewBarComponent",
          docgenInfo: OverviewBar_HTMLOverviewBarComponent.__docgenInfo,
          path: "src/components/xBars/OverviewBar.js"
      })
  },
  318: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__(15),
      __webpack_require__(24),
      __webpack_require__(21),
      __webpack_require__(5),
      __webpack_require__(26);
      var prop_types__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(1)
        , prop_types__WEBPACK_IMPORTED_MODULE_5___default = __webpack_require__.n(prop_types__WEBPACK_IMPORTED_MODULE_5__)
        , lodash_es__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(188)
        , lodash_es__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(320)
        , _CanvasComponent__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(159)
        , _store_connect__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(16)
        , _store_withPositionStore__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(37);
      function ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function _objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                  _defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }
      function _defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      class SequenceOverviewComponent extends _CanvasComponent__WEBPACK_IMPORTED_MODULE_8__.a {
          constructor(...args) {
              super(...args),
              _defineProperty(this, "draw", ()=>{
                  this.drawScene()
              }
              )
          }
          updateScrollPosition() {
              this._draw()
          }
          drawScene() {
              this.scene = {};
              var _this$props$position = this.props.position;
              this.scene.xViewPos = _this$props$position.xPos,
              this.scene.yViewPos = _this$props$position.yPos,
              this.scene.xScalingFactor = 1 / this.props.globalTileWidth * this.props.tileWidth,
              this.scene.yScalingFactor = 1 / this.props.globalTileHeight * this.props.tileHeight,
              this.drawCurrentViewpoint(),
              this.drawSequences()
          }
          drawSequences() {
              const _this$scene = this.scene
                , xViewPos = _this$scene.xViewPos
                , xScalingFactor = _this$scene.xScalingFactor
                , sequences = this.props.sequences.raw;
              let yPos = 0
                , i = 0;
              for (; i < sequences.length; i++) {
                  const sequence = sequences[i].sequence;
                  let xPos = 0
                    , j = Object(lodash_es__WEBPACK_IMPORTED_MODULE_6__.a)(Object(lodash_es__WEBPACK_IMPORTED_MODULE_7__.a)(xViewPos * xScalingFactor), 0, sequence.length - 1);
                  for (j = 0; j < sequence.length; j++) {
                      const el = sequence[j];
                      if (this.ctx.fillStyle(this.props.colorScheme.getColor(el)),
                      this.ctx.globalAlpha(.5),
                      this.ctx.fillRect(xPos, yPos, this.props.tileWidth, this.props.tileHeight),
                      xPos += this.props.tileWidth,
                      xPos > this.props.globalWidth)
                          break
                  }
                  if (yPos += this.props.tileHeight,
                  yPos > this.props.height)
                      break
              }
          }
          drawCurrentViewpoint() {
              const _this$scene2 = this.scene
                , xViewPos = _this$scene2.xViewPos
                , xScalingFactor = _this$scene2.xScalingFactor
                , yViewPos = _this$scene2.yViewPos
                , yScalingFactor = _this$scene2.yScalingFactor;
              this.ctx.globalAlpha(.8),
              this.ctx.fillRect(xViewPos * xScalingFactor, yViewPos * yScalingFactor, this.props.globalWidth * xScalingFactor, this.props.globalHeight * yScalingFactor)
          }
          render() {
              return super.render()
          }
      }
      SequenceOverviewComponent.displayName = "SequenceOverviewComponent",
      SequenceOverviewComponent.defaultProps = _objectSpread({}, _CanvasComponent__WEBPACK_IMPORTED_MODULE_8__.a.defaultProps, {
          height: 50,
          tileWidth: 5,
          tileHeight: 5
      }),
      SequenceOverviewComponent.propTypes = _objectSpread({}, _CanvasComponent__WEBPACK_IMPORTED_MODULE_8__.a.propTypes, {
          height: prop_types__WEBPACK_IMPORTED_MODULE_5___default.a.number,
          tileWidth: prop_types__WEBPACK_IMPORTED_MODULE_5___default.a.number,
          tileHeight: prop_types__WEBPACK_IMPORTED_MODULE_5___default.a.number
      });
      const SOC = Object(_store_withPositionStore__WEBPACK_IMPORTED_MODULE_10__.a)(SequenceOverviewComponent);
      __webpack_exports__.a = Object(_store_connect__WEBPACK_IMPORTED_MODULE_9__.a)(state=>({
          sequences: state.sequences,
          width: state.props.width,
          globalWidth: state.props.width,
          globalHeight: state.props.height,
          globalTileWidth: state.props.tileWidth,
          globalTileHeight: state.props.tileHeight,
          colorScheme: state.props.colorScheme
      }))(SOC),
      SequenceOverviewComponent.__docgenInfo = {
          description: "",
          methods: [{
              name: "draw",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "updateScrollPosition",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "drawScene",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "drawSequences",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "drawCurrentViewpoint",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }],
          displayName: "SequenceOverviewComponent",
          props: {
              height: {
                  defaultValue: {
                      value: "50",
                      computed: !1
                  },
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Height of the SequenceOverview (in pixels), e.g. `50`"
              },
              tileWidth: {
                  defaultValue: {
                      value: "5",
                      computed: !1
                  },
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Width of a tile in the OverviewBar, e.g. `5`"
              },
              tileHeight: {
                  defaultValue: {
                      value: "5",
                      computed: !1
                  },
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Height of a tile in the OverviewBar, e.g. `5`"
              }
          },
          composes: ["./CanvasComponent"]
      },
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/Canvas/SequenceOverview.js"] = {
          name: "SequenceOverviewComponent",
          docgenInfo: SequenceOverviewComponent.__docgenInfo,
          path: "src/components/Canvas/SequenceOverview.js"
      })
  },
  37: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__(15),
      __webpack_require__(21),
      __webpack_require__(5);
      var react__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(0)
        , react__WEBPACK_IMPORTED_MODULE_3___default = __webpack_require__.n(react__WEBPACK_IMPORTED_MODULE_3__)
        , prop_types__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(1)
        , prop_types__WEBPACK_IMPORTED_MODULE_4___default = __webpack_require__.n(prop_types__WEBPACK_IMPORTED_MODULE_4__)
        , lodash_es__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(319)
        , lodash_es__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(914)
        , _assert__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(117);
      function ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function _defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      __webpack_exports__.a = function withPositionConsumer(Component, {withX: withX=!1, withY: withY=!1}={}) {
          class MSAPositionConsumer extends react__WEBPACK_IMPORTED_MODULE_3__.PureComponent {
              constructor(props) {
                  super(props),
                  _defineProperty(this, "updateFromPositionStore", ()=>{
                      Object(_assert__WEBPACK_IMPORTED_MODULE_7__.a)(this.context && this.context.positionMSAStore, "MSA PositionStore needs to be injected");
                      const state = this.context.positionMSAStore.getState();
                      this.position = this.position || {};
                      const newPosition = Object(lodash_es__WEBPACK_IMPORTED_MODULE_5__.a)(state, ["currentViewSequence", "currentViewSequencePosition", "xPosOffset", "yPosOffset"]);
                      state.position && (newPosition.xPos = state.position.xPos,
                      newPosition.yPos = state.position.yPos),
                      this.setState({
                          highlight: state.highlight,
                          hasBeenInitialized: !0
                      }),
                      this.el.current && this.shouldRerender(newPosition) ? (this.position = newPosition,
                      this.setState({
                          position: this.position
                      })) : (Object(lodash_es__WEBPACK_IMPORTED_MODULE_6__.a)(newPosition, (v,k)=>{
                          this.position[k] = v
                      }
                      ),
                      this.el.current && this.updateScrollPosition())
                  }
                  ),
                  _defineProperty(this, "shouldRerender", newPosition=>{
                      const it = this.el.current;
                      if (void 0 !== it.shouldRerender)
                          return it.shouldRerender(newPosition);
                      const cacheElements = it.props.cacheElements;
                      return !!(withY && Math.abs(newPosition.currentViewSequence - this.position.lastCurrentViewSequence) >= cacheElements) || !!(withX && Math.abs(newPosition.currentViewSequencePosition - this.position.lastCurrentViewSequencePosition) >= cacheElements)
                  }
                  ),
                  _defineProperty(this, "updateScrollPosition", ()=>{
                      const it = this.el.current;
                      if (it && void 0 !== it.updateScrollPosition)
                          it.updateScrollPosition();
                      else if (it && it.el && it.el.current) {
                          if (withX) {
                              const tileWidth = it.props.tileWidth;
                              let offsetX = -this.position.xPosOffset;
                              offsetX += (this.position.lastCurrentViewSequencePosition - this.position.lastStartXTile) * tileWidth,
                              this.position.currentViewSequencePosition !== this.position.lastCurrentViewSequencePosition && (offsetX += (this.position.currentViewSequencePosition - this.position.lastCurrentViewSequencePosition) * tileWidth),
                              it.el.current.scrollLeft = offsetX
                          }
                          if (withY) {
                              const tileHeight = it.props.tileHeight;
                              let offsetY = -this.position.yPosOffset;
                              offsetY += (this.position.lastCurrentViewSequence - this.position.lastStartYTile) * tileHeight,
                              this.position.currentViewSequence !== this.position.lastCurrentViewSequence && (offsetY += (this.position.currentViewSequence - this.position.lastCurrentViewSequence) * tileHeight),
                              it.el.current.scrollTop = offsetY
                          }
                      }
                  }
                  ),
                  _defineProperty(this, "dispatch", payload=>{
                      this.context.positionMSAStore.dispatch(payload)
                  }
                  ),
                  this.el = react__WEBPACK_IMPORTED_MODULE_3___default.a.createRef(),
                  this.state = {
                      highlight: null,
                      hasBeenInitialized: !1
                  }
              }
              componentDidMount() {
                  this.unsubscribe = this.context.positionMSAStore.subscribe(this.updateFromPositionStore),
                  this.updateScrollPosition(!0),
                  this.updateFromPositionStore()
              }
              componentDidUpdate() {
                  this.updateScrollPosition()
              }
              componentWillUnmount() {
                  this.unsubscribe()
              }
              render() {
                  return this.state.hasBeenInitialized ? react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(Component, function _objectSpread(target) {
                      for (var source, i = 1; i < arguments.length; i++)
                          source = null != arguments[i] ? arguments[i] : {},
                          i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                              _defineProperty(target, key, source[key])
                          }
                          )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                              Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
                          }
                          ));
                      return target
                  }({
                      ref: this.el,
                      position: this.position,
                      positionDispatch: this.dispatch,
                      highlight: this.state.highlight
                  }, this.props)) : null
              }
          }
          return MSAPositionConsumer.displayName = "MSAPositionConsumer",
          MSAPositionConsumer.displayName = "withPosition(".concat(Component.displayName || Component.name, ")"),
          MSAPositionConsumer.contextTypes = {
              positionMSAStore: prop_types__WEBPACK_IMPORTED_MODULE_4___default.a.object
          },
          MSAPositionConsumer
      }
  },
  47: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.d(__webpack_exports__, "b", (function() {
          return movePosition
      }
      )),
      __webpack_require__.d(__webpack_exports__, "a", (function() {
          return actions
      }
      )),
      __webpack_require__.d(__webpack_exports__, "c", (function() {
          return positionReducer
      }
      ));
      __webpack_require__(15),
      __webpack_require__(21),
      __webpack_require__(5);
      var _actions__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(8)
        , lodash_es__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(188)
        , lodash_es__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(319)
        , lodash_es__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(320)
        , _assert__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(117);
      function ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function _objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                  _defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }
      function _defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      const updateMainStore = Object(_actions__WEBPACK_IMPORTED_MODULE_3__.a)("MAINSTORE_UPDATE")
        , movePosition = Object(_actions__WEBPACK_IMPORTED_MODULE_3__.a)("POSITION_MOVE")
        , updatePosition = Object(_actions__WEBPACK_IMPORTED_MODULE_3__.a)("POSITION_UPDATE")
        , updatePositionByResidue = Object(_actions__WEBPACK_IMPORTED_MODULE_3__.a)("POSITION_BY_RESIDUE_UPDATE")
        , highlightRegion = Object(_actions__WEBPACK_IMPORTED_MODULE_3__.a)("HIGHLIGHT_REGION")
        , removeHighlightRegion = Object(_actions__WEBPACK_IMPORTED_MODULE_3__.a)("REMOVE_HIGHLIGHT_REGION")
        , actions = {
          updateMainStore: updateMainStore,
          updatePosition: updatePosition,
          movePosition: movePosition,
          updatePositionByResidue: updatePositionByResidue,
          highlightRegion: highlightRegion,
          removeHighlightRegion: removeHighlightRegion
      };
      function commonPositionReducer(prevState, pos) {
          const maxWidth = prevState.sequences.maxLength * prevState.props.tileWidth - prevState.props.width;
          pos.xPos = Object(lodash_es__WEBPACK_IMPORTED_MODULE_4__.a)(pos.xPos, 0, maxWidth);
          const maxHeight = prevState.sequences.raw.length * prevState.props.tileHeight - prevState.props.height;
          return pos.yPos = Object(lodash_es__WEBPACK_IMPORTED_MODULE_4__.a)(pos.yPos, 0, maxHeight),
          _objectSpread({}, prevState, {
              position: pos
          })
      }
      function positionReducer(oldState={
          position: {
              xPos: 0,
              yPos: 0
          },
          highlight: null
      }, action) {
          let state = oldState
            , position = oldState.position
            , highlight = oldState.highlight;
          switch (action.type) {
          case updateMainStore.key:
              state = _objectSpread({}, Object(lodash_es__WEBPACK_IMPORTED_MODULE_5__.a)(state, ["props", "sequenceStats", "sequences"]), {}, action.payload);
              break;
          case updatePosition.key:
          case movePosition.key:
          case updatePositionByResidue.key:
              position = ((prevState={
                  position: {
                      xPos: 0,
                      yPos: 0
                  }
              },action)=>{
                  const pos = prevState.position;
                  switch (action.type) {
                  case movePosition.key:
                      Object(_assert__WEBPACK_IMPORTED_MODULE_7__.a)(void 0 !== action.payload.xMovement || void 0 !== action.payload.yMovement, "must contain at least xMovement or yMovement");
                      const movePayload = _objectSpread({}, pos);
                      return movePayload.xPos += action.payload.xMovement || 0,
                      movePayload.yPos += action.payload.yMovement || 0,
                      commonPositionReducer(prevState, movePayload);
                  case updatePosition.key:
                      return Object(_assert__WEBPACK_IMPORTED_MODULE_7__.a)(void 0 !== action.payload.xPos || void 0 !== action.payload.yPos, "must contain at least xPos or yPos"),
                      commonPositionReducer(prevState, {
                          xPos: isNaN(action.payload.xPos) ? pos.xPos : action.payload.xPos,
                          yPos: isNaN(action.payload.yPos) ? pos.yPos : action.payload.yPos
                      });
                  case updatePositionByResidue.key:
                      return Object(_assert__WEBPACK_IMPORTED_MODULE_7__.a)(void 0 !== action.payload.aaPos, "must contain at least aaPos"),
                      commonPositionReducer(prevState, {
                          xPos: prevState.props.tileWidth * (action.payload.aaPos - 1),
                          yPos: pos.yPos
                      });
                  default:
                      return prevState
                  }
              }
              )(state, action).position;
              break;
          case highlightRegion.key:
              highlight = action.payload;
              break;
          case removeHighlightRegion.key:
              highlight = null;
              break;
          default:
              return state
          }
          return _objectSpread({}, state, {}, {
              xPosOffset: -position.xPos % state.props.tileWidth,
              yPosOffset: -position.yPos % state.props.tileWidth,
              currentViewSequence: Object(lodash_es__WEBPACK_IMPORTED_MODULE_4__.a)(Object(lodash_es__WEBPACK_IMPORTED_MODULE_6__.a)(position.yPos / state.props.tileHeight), 0, state.sequences.length - 1),
              currentViewSequencePosition: Object(lodash_es__WEBPACK_IMPORTED_MODULE_4__.a)(Object(lodash_es__WEBPACK_IMPORTED_MODULE_6__.a)(position.xPos / state.props.tileWidth), 0, state.sequences.maxLength),
              position: position
          }, {
              highlight: highlight
          })
      }
  },
  470: function(module, exports, __webpack_require__) {
      __webpack_require__(471),
      __webpack_require__(616),
      module.exports = __webpack_require__(617)
  },
  535: function(module, exports) {},
  57: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.d(__webpack_exports__, "d", (function() {
          return SequencePropType
      }
      )),
      __webpack_require__.d(__webpack_exports__, "a", (function() {
          return ColorSchemePropType
      }
      )),
      __webpack_require__.d(__webpack_exports__, "c", (function() {
          return PositionPropType
      }
      )),
      __webpack_require__.d(__webpack_exports__, "b", (function() {
          return MSAPropTypes
      }
      )),
      __webpack_require__.d(__webpack_exports__, "e", (function() {
          return msaDefaultProps
      }
      ));
      var prop_types__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(1)
        , prop_types__WEBPACK_IMPORTED_MODULE_0___default = __webpack_require__.n(prop_types__WEBPACK_IMPORTED_MODULE_0__)
        , _utils_ColorScheme__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(112);
      const SequencePropType = prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.shape({
          name: prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.string,
          sequence: prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.string
      })
        , ColorSchemePropType = prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.oneOfType([prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.oneOf(["aliphatic", "aromatic", "charged", "negative", "polar", "positive", "serine_threonine", "buried_index", "clustal", "clustal2", "cinema", "helix_propensity", "hydro", "lesk", "mae", "nucleotide", "purine_pyrimidine", "strand_propensity", "taylor", "turn_propensity", "zappo", "conservation"]), prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.instanceOf(_utils_ColorScheme__WEBPACK_IMPORTED_MODULE_1__.a), function(props, propName, componentName) {
          if (!Object(_utils_ColorScheme__WEBPACK_IMPORTED_MODULE_1__.b)(props[propName]))
              return new Error("Invalid prop `" + propName + "` supplied to `" + componentName + "`. Validation failed.")
      }
      ])
        , PositionPropType = prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.shape({
          xPos: prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.number,
          yPos: prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.number
      })
        , MSAPropTypes = {
          width: prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.number,
          height: prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.number,
          tileWidth: prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.number,
          tileHeight: prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.number,
          position: PositionPropType,
          highlight: prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.oneOfType([prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.object, prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.array]),
          colorScheme: ColorSchemePropType,
          calculateConservation: prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.bool,
          overlayConservation: prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.bool,
          sampleSizeConservation: prop_types__WEBPACK_IMPORTED_MODULE_0___default.a.number
      }
        , msaDefaultProps = {
          width: 800,
          height: 600,
          tileWidth: 20,
          tileHeight: 20,
          colorScheme: "clustal",
          calculateConservation: !1,
          overlayConservation: !1,
          sampleSizeConservation: void 0
      }
  },
  617: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.r(__webpack_exports__),
      function(module) {
          var _storybook_react__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(31)
            , _storybook_addon_options__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(462);
          Object(_storybook_react__WEBPACK_IMPORTED_MODULE_0__.addDecorator)(Object(_storybook_addon_options__WEBPACK_IMPORTED_MODULE_1__.withOptions)({
              name: "React MSA Viewer",
              url: "#",
              goFullScreen: !1,
              showStoriesPanel: !0,
              showAddonPanel: !0,
              showSearchBox: !1,
              addonPanelInRight: !1,
              sortStoriesByKind: !1,
              hierarchySeparator: null,
              hierarchyRootSeparator: null,
              sidebarAnimations: !0,
              selectedAddonPanel: void 0,
              enableShortcuts: !1
          })),
          Object(_storybook_react__WEBPACK_IMPORTED_MODULE_0__.configure)((function loadStories() {
              __webpack_require__(851)
          }
          ), module)
      }
      .call(this, __webpack_require__(49)(module))
  },
  76: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      var reselect__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(185)
        , _shallowEqual__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(93);
      const shallowSelect = Object(reselect__WEBPACK_IMPORTED_MODULE_0__.b)(reselect__WEBPACK_IMPORTED_MODULE_0__.c, _shallowEqual__WEBPACK_IMPORTED_MODULE_1__.a);
      __webpack_exports__.a = shallowSelect
  },
  77: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.d(__webpack_exports__, "a", (function() {
          return autoBind
      }
      ));
      var _assert__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(117);
      function autoBind(instance) {
          Object(_assert__WEBPACK_IMPORTED_MODULE_0__.a)(1 < arguments.length, "Must provide methods for binding");
          for (let i = 1; i < arguments.length; i++) {
              const k = arguments[i];
              instance[k] = instance[k].bind(instance)
          }
      }
  },
  8: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.d(__webpack_exports__, "a", (function() {
          return createAction
      }
      )),
      __webpack_require__.d(__webpack_exports__, "c", (function() {
          return updateProps
      }
      )),
      __webpack_require__.d(__webpack_exports__, "d", (function() {
          return updateSequences
      }
      ));
      __webpack_require__(52),
      __webpack_require__(5);
      function createAction(type, ...argNames) {
          const actionCreator = function(...args) {
              let payload;
              return 0 === argNames.length ? payload = args[0] : (payload = {},
              argNames.forEach((arg,index)=>{
                  payload[argNames[index]] = args[index]
              }
              )),
              {
                  type: type,
                  payload: payload
              }
          };
          return actionCreator.toString = ()=>type.toString(),
          actionCreator.key = actionCreator.toString(),
          actionCreator
      }
      const updateProps = createAction("PROPS_UPDATE")
        , updateProp = createAction("PROP_UPDATE", "key", "value")
        , updateSequences = createAction("SEQUENCES_UPDATE")
        , updateConservation = createAction("CONSERVSTION_UPDATE")
        , actions = {
          updateProp: updateProp,
          updateProps: updateProps,
          updateSequences: updateSequences,
          updateConservation: updateConservation
      };
      __webpack_exports__.b = actions
  },
  851: function(module, exports, __webpack_require__) {
      __webpack_require__(852),
      __webpack_require__(895),
      __webpack_require__(896),
      __webpack_require__(897),
      __webpack_require__(898),
      __webpack_require__(911),
      __webpack_require__(912)
  },
  852: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.r(__webpack_exports__),
      function(module) {
          var react__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(0)
            , react__WEBPACK_IMPORTED_MODULE_0___default = __webpack_require__.n(react__WEBPACK_IMPORTED_MODULE_0__)
            , _storybook_react__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(31)
            , _lib__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(10)
            , lodash_es__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(317)
            , _storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(9)
            , withSourceLoader = __webpack_require__(38).withSource
            , __STORY__ = (__webpack_require__(38).addSource,
          '/**\n * Copyright 2018, Plotly, Inc.\n * All rights reserved.\n *\n * This source code is licensed under the MIT license found in the\n * LICENSE file in the root directory of this source tree.\n */\n\nimport React from "react";\nimport { storiesOf } from "@storybook/react";\nimport { MSAViewer } from "../lib";\nimport { times } from "lodash-es";\nimport { number, withKnobs } from "@storybook/addon-knobs";\n\nconst sequences = [\n  {\n    name: "seq.1",\n    sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED",\n  },\n  {\n    name: "seq.2",\n    sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n  {\n    name: "seq.3",\n    sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n  {\n    name: "seq.4",\n    sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n  {\n    name: "seq.5",\n    sequence: "MEEPQSD--IEL-PLSEETFSDLWWPLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n  {\n    name: "seq.6",\n    sequence: "MEEPQEDLSSSL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n  {\n    name: "seq.7",\n    sequence: "MEEPQ---SISE-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE---LSENVAGWLEDP",\n  },\n];\n\nstoriesOf("Basic", module)\n  .add("Standard rendering", function () {\n    const options = {\n      sequences,\n    };\n    return <MSAViewer {...options} />;\n  })\n  .addDecorator(withKnobs)\n  .add("Big viewpoint", function () {\n    const options = {\n      sequences: [],\n      height: number("height", 500),\n      width: number("width", 500),\n      tileHeight: number("tileHeight", 20),\n      tileWidth: number("tileWidth", 20),\n      colorScheme: "clustal",\n    };\n    const sequence =\n      "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED";\n    times(100, () => {\n      const mutation_pos = Math.round(Math.random() * sequence.length);\n      options.sequences.push({\n        name: `seq_${options.sequences.length}`,\n        sequence:\n          sequence.substring(0, mutation_pos) +\n          "-" +\n          sequence.substring(mutation_pos + 1),\n      });\n    });\n    return <MSAViewer {...options} />;\n  });\n')
            , __ADDS_MAP__ = {
              "basic--big-viewpoint": {
                  startLoc: {
                      col: 7,
                      line: 54
                  },
                  endLoc: {
                      col: 3,
                      line: 76
                  },
                  startBody: {
                      col: 24,
                      line: 54
                  },
                  endBody: {
                      col: 3,
                      line: 76
                  }
              },
              "basic--standard-rendering": {
                  startLoc: {
                      col: 7,
                      line: 47
                  },
                  endLoc: {
                      col: 3,
                      line: 52
                  },
                  startBody: {
                      col: 29,
                      line: 47
                  },
                  endBody: {
                      col: 3,
                      line: 52
                  }
              }
          };
          const sequences = [{
              name: "seq.1",
              sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED"
          }, {
              name: "seq.2",
              sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
          }, {
              name: "seq.3",
              sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
          }, {
              name: "seq.4",
              sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
          }, {
              name: "seq.5",
              sequence: "MEEPQSD--IEL-PLSEETFSDLWWPLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
          }, {
              name: "seq.6",
              sequence: "MEEPQEDLSSSL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
          }, {
              name: "seq.7",
              sequence: "MEEPQ---SISE-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE---LSENVAGWLEDP"
          }];
          Object(_storybook_react__WEBPACK_IMPORTED_MODULE_1__.storiesOf)("Basic", module).addParameters({
              storySource: {
                  source: __STORY__,
                  locationsMap: __ADDS_MAP__
              }
          }).addDecorator(withSourceLoader(__STORY__, __ADDS_MAP__, "/basic.js", [], {}, "/opt/build/repo/src/stories", {})).add("Standard rendering", (function() {
              return react__WEBPACK_IMPORTED_MODULE_0___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_2__.a, {
                  sequences: sequences
              })
          }
          )).addDecorator(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_4__.withKnobs).add("Big viewpoint", (function() {
              const options = {
                  sequences: [],
                  height: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_4__.number)("height", 500),
                  width: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_4__.number)("width", 500),
                  tileHeight: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_4__.number)("tileHeight", 20),
                  tileWidth: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_4__.number)("tileWidth", 20),
                  colorScheme: "clustal"
              }
                , sequence = "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED";
              return Object(lodash_es__WEBPACK_IMPORTED_MODULE_3__.a)(100, ()=>{
                  const mutation_pos = Math.round(Math.random() * sequence.length);
                  options.sequences.push({
                      name: "seq_".concat(options.sequences.length),
                      sequence: sequence.substring(0, mutation_pos) + "-" + sequence.substring(mutation_pos + 1)
                  })
              }
              ),
              react__WEBPACK_IMPORTED_MODULE_0___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_2__.a, options)
          }
          ))
      }
      .call(this, __webpack_require__(49)(module))
  },
  895: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.r(__webpack_exports__),
      function(module) {
          __webpack_require__(45),
          __webpack_require__(33),
          __webpack_require__(103),
          __webpack_require__(24),
          __webpack_require__(178),
          __webpack_require__(52),
          __webpack_require__(26);
          var react__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(0)
            , react__WEBPACK_IMPORTED_MODULE_7___default = __webpack_require__.n(react__WEBPACK_IMPORTED_MODULE_7__)
            , _storybook_react__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(31)
            , _lib__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(10)
            , lodash_es__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(317)
            , _storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_11__ = __webpack_require__(9);
          function _createForOfIteratorHelper(o) {
              if ("undefined" == typeof Symbol || null == o[Symbol.iterator]) {
                  if (Array.isArray(o) || (o = function _unsupportedIterableToArray(o, minLen) {
                      if (!o)
                          return;
                      if ("string" == typeof o)
                          return _arrayLikeToArray(o, minLen);
                      var n = Object.prototype.toString.call(o).slice(8, -1);
                      "Object" === n && o.constructor && (n = o.constructor.name);
                      if ("Map" === n || "Set" === n)
                          return Array.from(n);
                      if ("Arguments" === n || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n))
                          return _arrayLikeToArray(o, minLen)
                  }(o))) {
                      var i = 0
                        , F = function() {};
                      return {
                          s: F,
                          n: function n() {
                              return i >= o.length ? {
                                  done: !0
                              } : {
                                  done: !1,
                                  value: o[i++]
                              }
                          },
                          e: function e(_e) {
                              throw _e
                          },
                          f: F
                      }
                  }
                  throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.")
              }
              var it, err, normalCompletion = !0, didErr = !1;
              return {
                  s: function s() {
                      it = o[Symbol.iterator]()
                  },
                  n: function n() {
                      var step = it.next();
                      return normalCompletion = step.done,
                      step
                  },
                  e: function e(_e2) {
                      didErr = !0,
                      err = _e2
                  },
                  f: function f() {
                      try {
                          normalCompletion || null == it.return || it.return()
                      } finally {
                          if (didErr)
                              throw err
                      }
                  }
              }
          }
          function _arrayLikeToArray(arr, len) {
              (null == len || len > arr.length) && (len = arr.length);
              for (var i = 0, arr2 = Array(len); i < len; i++)
                  arr2[i] = arr[i];
              return arr2
          }
          var withSourceLoader = __webpack_require__(38).withSource
            , __STORY__ = (__webpack_require__(38).addSource,
          'import React, { Component } from "react";\nimport { storiesOf } from "@storybook/react";\nimport { MSAViewer } from "../lib";\nimport { times } from "lodash-es";\nimport { number } from "@storybook/addon-knobs";\n\nconst alphabet = "ACDEFGHIKLMNPQRSTVWY";\nconst getRandomBase = () =>\n  alphabet[Math.floor(Math.random() * alphabet.length)];\nconst generateSequence = (length) => {\n  let seq = "";\n  for (let i = 0; i < length; i++) seq += getRandomBase();\n  return seq;\n};\nconst seqLengths = [100, 1000, 10000];\nconst nSeqs = [100, 1000, 10000];\n\nfor (let seqLength of seqLengths) {\n  for (let nSeq of nSeqs) {\n    storiesOf("Performance Test", module).add(\n      `${nSeq} sequences of ${seqLength} residues`,\n      function () {\n        const options = {\n          sequences: [],\n          height: number("height", 500),\n          width: number("width", 500),\n          tileHeight: number("tileHeight", 20),\n          tileWidth: number("tileWidth", 20),\n          colorScheme: "clustal",\n        };\n        let time = Date.now();\n        const sequence = generateSequence(seqLength);\n        times(nSeq, () => {\n          const mutation_pos = Math.round(\n            Math.random() * (sequence.length - 1)\n          );\n          options.sequences.push({\n            name: `seq_${options.sequences.length}`,\n            sequence:\n              sequence.substring(0, mutation_pos) +\n              "-" +\n              sequence.substring(mutation_pos + 1),\n          });\n        });\n        time = Date.now() - time;\n        const resetView = Date.now();\n\n        class ViewerWithPerformance extends Component {\n          state = { resetView: null };\n          componentDidMount() {\n            this.setState({ resetView: Date.now() - resetView });\n          }\n          render() {\n            return (\n              <>\n                <h5>Performance time in ms</h5>\n                <code>\n                  generateSequnces => {time} || resetView =>{" "}\n                  {this.state.resetView}\n                </code>\n                <MSAViewer {...options} />\n              </>\n            );\n          }\n        }\n        return <ViewerWithPerformance />;\n      }\n    );\n  }\n}\n')
            , __ADDS_MAP__ = {}
            , __MODULE_DEPENDENCIES__ = []
            , __LOCAL_DEPENDENCIES__ = {}
            , __IDS_TO_FRAMEWORKS__ = {};
          const generateSequence = length=>{
              let seq = "";
              for (let i = 0; i < length; i++)
                  seq += "ACDEFGHIKLMNPQRSTVWY"[Math.floor(Math.random() * "ACDEFGHIKLMNPQRSTVWY".length)];
              return seq
          }
            , seqLengths = [100, 1e3, 1e4]
            , nSeqs = [100, 1e3, 1e4];
          for (var _ref = react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("h5", null, "Performance time in ms"), _i = 0, _seqLengths = seqLengths; _i < _seqLengths.length; _i++) {
              let seqLength = _seqLengths[_i];
              var _step, _iterator = _createForOfIteratorHelper(nSeqs);
              try {
                  for (_iterator.s(); !(_step = _iterator.n()).done; ) {
                      let nSeq = _step.value;
                      Object(_storybook_react__WEBPACK_IMPORTED_MODULE_8__.storiesOf)("Performance Test", module).addParameters({
                          storySource: {
                              source: __STORY__,
                              locationsMap: __ADDS_MAP__
                          }
                      }).addDecorator(withSourceLoader(__STORY__, __ADDS_MAP__, "/performanceTest.js", __MODULE_DEPENDENCIES__, __LOCAL_DEPENDENCIES__, "/opt/build/repo/src/stories", __IDS_TO_FRAMEWORKS__)).add("".concat(nSeq, " sequences of ").concat(seqLength, " residues"), (function() {
                          const options = {
                              sequences: [],
                              height: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_11__.number)("height", 500),
                              width: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_11__.number)("width", 500),
                              tileHeight: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_11__.number)("tileHeight", 20),
                              tileWidth: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_11__.number)("tileWidth", 20),
                              colorScheme: "clustal"
                          };
                          let time = Date.now();
                          const sequence = generateSequence(seqLength);
                          Object(lodash_es__WEBPACK_IMPORTED_MODULE_10__.a)(nSeq, ()=>{
                              const mutation_pos = Math.round(Math.random() * (sequence.length - 1));
                              options.sequences.push({
                                  name: "seq_".concat(options.sequences.length),
                                  sequence: sequence.substring(0, mutation_pos) + "-" + sequence.substring(mutation_pos + 1)
                              })
                          }
                          ),
                          time = Date.now() - time;
                          const resetView = Date.now();
                          class ViewerWithPerformance extends react__WEBPACK_IMPORTED_MODULE_7__.Component {
                              constructor(...args) {
                                  var obj, key, value;
                                  super(...args),
                                  value = {
                                      resetView: null
                                  },
                                  (key = "state")in (obj = this) ? Object.defineProperty(obj, key, {
                                      value: value,
                                      enumerable: !0,
                                      configurable: !0,
                                      writable: !0
                                  }) : obj[key] = value
                              }
                              componentDidMount() {
                                  this.setState({
                                      resetView: Date.now() - resetView
                                  })
                              }
                              render() {
                                  return react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(react__WEBPACK_IMPORTED_MODULE_7___default.a.Fragment, null, _ref, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("code", null, "generateSequnces => ", time, " || resetView =>", " ", this.state.resetView), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_9__.a, options))
                              }
                          }
                          return ViewerWithPerformance.displayName = "ViewerWithPerformance",
                          react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(ViewerWithPerformance, null)
                      }
                      ))
                  }
              } catch (err) {
                  _iterator.e(err)
              } finally {
                  _iterator.f()
              }
          }
      }
      .call(this, __webpack_require__(49)(module))
  },
  896: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.r(__webpack_exports__),
      function(module) {
          __webpack_require__(24),
          __webpack_require__(35),
          __webpack_require__(7),
          __webpack_require__(26);
          var react__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(0)
            , react__WEBPACK_IMPORTED_MODULE_4___default = __webpack_require__.n(react__WEBPACK_IMPORTED_MODULE_4__)
            , _storybook_react__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(31)
            , _lib__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(10)
            , _src_colorschemes__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(187)
            , _storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(9)
            , lodash_es__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(916);
          function _extends() {
              return (_extends = Object.assign || function(target) {
                  for (var source, i = 1; i < arguments.length; i++)
                      for (var key in source = arguments[i])
                          Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
                  return target
              }
              ).apply(this, arguments)
          }
          var withSourceLoader = __webpack_require__(38).withSource
            , __STORY__ = (__webpack_require__(38).addSource,
          '/**\n * Copyright 2018, Plotly, Inc.\n * All rights reserved.\n *\n * This source code is licensed under the MIT license found in the\n * LICENSE file in the root directory of this source tree.\n */\n\nimport React, { Component } from "react";\nimport { storiesOf } from "@storybook/react";\nimport { MSAViewer } from "../lib";\nimport { staticSchemes, dynSchemes } from "../../src/colorschemes";\nimport {\n  select,\n  text,\n  boolean,\n  number,\n  withKnobs,\n} from "@storybook/addon-knobs";\nimport { zipObject } from "lodash-es";\n\nconst sequences = [\n  {\n    name: "seq.1",\n    sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED",\n  },\n  {\n    name: "seq.2",\n    sequence: "MEEPQSDLSIEL-PLSQETFTDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n  {\n    name: "seq.3",\n    sequence: "MEPIQSDLSIEL-PLSQETFWDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n  {\n    name: "seq.4",\n    sequence: "MIPEQSSLSIEL-PLSQETFLDLWKLYPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n];\n\n// Storybook 4 selects only accepts key/value objects\nfunction createObject(options) {\n  return zipObject(options, options);\n}\n\nstoriesOf("Customization", module)\n  .addDecorator(withKnobs)\n  .add("Colorschemes", function () {\n    // 2018: see https://github.com/wilzbach/msa-colorschemes for now\n    // 2020: additional color schemes added for UniProt amino acid properties\n    const colorSchemes = [\n      ...Object.keys(staticSchemes),\n      ...Object.keys(dynSchemes),\n    ];\n    const options = {\n      colorScheme: select("Colorscheme", createObject(colorSchemes), "zappo"),\n      calculateConservation: true,\n      overlayConservation: boolean("overlayConservation", false),\n      sampleSizeConservation: number("sampleSizeConservation", ""),\n      sequenceTextFont: "16px Monospace",\n      sequences,\n    };\n    let currentColor = null;\n    class MSAConservation extends Component {\n      componentDidMount() {\n        this.parentDiv.addEventListener("conservationProgress", console.log);\n        this.parentDiv.addEventListener("drawCompleted", () => {\n          const { name, map } = this.viewer.getColorMap();\n          if (name !== currentColor) {\n            console.log(name, map);\n            currentColor = name;\n          }\n        });\n      }\n      render() {\n        return (\n          <div ref={(ref) => (this.parentDiv = ref)}>\n            <MSAViewer {...options} ref={(ref) => (this.viewer = ref)} />\n          </div>\n        );\n      }\n    }\n    return <MSAConservation />;\n  })\n  .add("Custom ColorScheme", function () {\n    // see https://github.com/wilzbach/msa-colorschemes for now\n    const myColorMap = {\n      M: "blue",\n      E: "red",\n      T: "green",\n    };\n    class MyColorScheme {\n      getColor(element) {\n        return element in myColorMap ? myColorMap[element] : "grey";\n      }\n    }\n    const myColorScheme = new MyColorScheme();\n    const options = {\n      colorScheme: myColorScheme,\n      sequences,\n    };\n    return <MSAViewer {...options} />;\n  })\n  .add("Custom Labels", function () {\n    const fontSizes = [\n      "6px",\n      "8px",\n      "10px",\n      "12px",\n      "14px",\n      "16px",\n      "18px",\n      "20px",\n    ];\n    const fontSize = select("Font size", createObject(fontSizes), "14px");\n    const options = {\n      sequences,\n      labelComponent: ({ sequence }) => {\n        return (\n          <div style={{ height: 20, fontWeight: "bold", fontSize }}>\n            My: {sequence.name}\n          </div>\n        );\n      },\n    };\n    return <MSAViewer {...options} />;\n  })\n  .add("Custom Markers", function () {\n    const fontSizes = ["6px", "8px", "10px", "12px", "14px", "16px", "18px"];\n    const fontSize = select("Font size", createObject(fontSizes), "10px");\n    const options = {\n      sequences,\n      markerComponent: ({ index }) => {\n        return (\n          <div\n            style={{\n              width: 20,\n              display: "inline-block",\n              textAlign: "center",\n              fontSize: fontSize,\n              fontWeight: "bold",\n              color: "pink",\n            }}\n          >\n            {index}\n          </div>\n        );\n      },\n    };\n    return <MSAViewer {...options} />;\n  })\n  .add("Custom styling", function () {\n    const options = {\n      sequences,\n      labelStyle: {\n        outline: text("Label style (outline)", "1px solid black"),\n      },\n      markerStyle: {\n        outline: text("Marker style (outline)", "1px solid black"),\n      },\n      sequenceTextColor: text("Sequence color", "blue"),\n    };\n    return <MSAViewer {...options} />;\n  })\n  .add("Custom scollbars", function () {\n    const options = {\n      sequences,\n      sequenceScrollBarPositionX: select(\n        "ScrollBarPositionX",\n        createObject(["top", "bottom"]),\n        "top"\n      ),\n      sequenceScrollBarPositionY: select(\n        "ScrollBarPositionY",\n        createObject(["left", "right"]),\n        "left"\n      ),\n      sequenceOverflow: select(\n        "Overflow",\n        createObject(["scroll", "auto", "hidden"]),\n        "scroll"\n      ),\n    };\n    return <MSAViewer {...options} />;\n  });\n')
            , __ADDS_MAP__ = {
              "customization--custom-scollbars": {
                  startLoc: {
                      col: 7,
                      line: 165
                  },
                  endLoc: {
                      col: 3,
                      line: 185
                  },
                  startBody: {
                      col: 27,
                      line: 165
                  },
                  endBody: {
                      col: 3,
                      line: 185
                  }
              },
              "customization--custom-styling": {
                  startLoc: {
                      col: 7,
                      line: 152
                  },
                  endLoc: {
                      col: 3,
                      line: 164
                  },
                  startBody: {
                      col: 25,
                      line: 152
                  },
                  endBody: {
                      col: 3,
                      line: 164
                  }
              },
              "customization--custom-markers": {
                  startLoc: {
                      col: 7,
                      line: 128
                  },
                  endLoc: {
                      col: 3,
                      line: 151
                  },
                  startBody: {
                      col: 25,
                      line: 128
                  },
                  endBody: {
                      col: 3,
                      line: 151
                  }
              },
              "customization--custom-labels": {
                  startLoc: {
                      col: 7,
                      line: 104
                  },
                  endLoc: {
                      col: 3,
                      line: 127
                  },
                  startBody: {
                      col: 24,
                      line: 104
                  },
                  endBody: {
                      col: 3,
                      line: 127
                  }
              },
              "customization--custom-colorscheme": {
                  startLoc: {
                      col: 7,
                      line: 85
                  },
                  endLoc: {
                      col: 3,
                      line: 103
                  },
                  startBody: {
                      col: 29,
                      line: 85
                  },
                  endBody: {
                      col: 3,
                      line: 103
                  }
              },
              "customization--colorschemes": {
                  startLoc: {
                      col: 7,
                      line: 48
                  },
                  endLoc: {
                      col: 3,
                      line: 84
                  },
                  startBody: {
                      col: 23,
                      line: 48
                  },
                  endBody: {
                      col: 3,
                      line: 84
                  }
              }
          };
          const sequences = [{
              name: "seq.1",
              sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED"
          }, {
              name: "seq.2",
              sequence: "MEEPQSDLSIEL-PLSQETFTDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
          }, {
              name: "seq.3",
              sequence: "MEPIQSDLSIEL-PLSQETFWDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
          }, {
              name: "seq.4",
              sequence: "MIPEQSSLSIEL-PLSQETFLDLWKLYPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
          }];
          function createObject(options) {
              return Object(lodash_es__WEBPACK_IMPORTED_MODULE_9__.a)(options, options)
          }
          Object(_storybook_react__WEBPACK_IMPORTED_MODULE_5__.storiesOf)("Customization", module).addParameters({
              storySource: {
                  source: __STORY__,
                  locationsMap: __ADDS_MAP__
              }
          }).addDecorator(withSourceLoader(__STORY__, __ADDS_MAP__, "/customization.js", [], {}, "/opt/build/repo/src/stories", {})).addDecorator(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__.withKnobs).add("Colorschemes", (function() {
              const colorSchemes = [...Object.keys(_src_colorschemes__WEBPACK_IMPORTED_MODULE_7__.c), ...Object.keys(_src_colorschemes__WEBPACK_IMPORTED_MODULE_7__.b)]
                , options = {
                  colorScheme: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__.select)("Colorscheme", createObject(colorSchemes), "zappo"),
                  calculateConservation: !0,
                  overlayConservation: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__.boolean)("overlayConservation", !1),
                  sampleSizeConservation: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__.number)("sampleSizeConservation", ""),
                  sequenceTextFont: "16px Monospace",
                  sequences: sequences
              };
              let currentColor = null;
              class MSAConservation extends react__WEBPACK_IMPORTED_MODULE_4__.Component {
                  componentDidMount() {
                      this.parentDiv.addEventListener("conservationProgress", console.log),
                      this.parentDiv.addEventListener("drawCompleted", ()=>{
                          const _this$viewer$getColor = this.viewer.getColorMap()
                            , name = _this$viewer$getColor.name
                            , map = _this$viewer$getColor.map;
                          name !== currentColor && (console.log(name, map),
                          currentColor = name)
                      }
                      )
                  }
                  render() {
                      return react__WEBPACK_IMPORTED_MODULE_4___default.a.createElement("div", {
                          ref: _ref2=>this.parentDiv = _ref2
                      }, react__WEBPACK_IMPORTED_MODULE_4___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_6__.a, _extends({}, options, {
                          ref: _ref=>this.viewer = _ref
                      })))
                  }
              }
              return MSAConservation.displayName = "MSAConservation",
              react__WEBPACK_IMPORTED_MODULE_4___default.a.createElement(MSAConservation, null)
          }
          )).add("Custom ColorScheme", (function() {
              const myColorMap = {
                  M: "blue",
                  E: "red",
                  T: "green"
              }
                , myColorScheme = new class {
                  getColor(element) {
                      return element in myColorMap ? myColorMap[element] : "grey"
                  }
              }
              ;
              return react__WEBPACK_IMPORTED_MODULE_4___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_6__.a, {
                  colorScheme: myColorScheme,
                  sequences: sequences
              })
          }
          )).add("Custom Labels", (function() {
              const fontSize = Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__.select)("Font size", createObject(["6px", "8px", "10px", "12px", "14px", "16px", "18px", "20px"]), "14px");
              return react__WEBPACK_IMPORTED_MODULE_4___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_6__.a, {
                  sequences: sequences,
                  labelComponent: ({sequence: sequence})=>react__WEBPACK_IMPORTED_MODULE_4___default.a.createElement("div", {
                      style: {
                          height: 20,
                          fontWeight: "bold",
                          fontSize: fontSize
                      }
                  }, "My: ", sequence.name)
              })
          }
          )).add("Custom Markers", (function() {
              const fontSize = Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__.select)("Font size", createObject(["6px", "8px", "10px", "12px", "14px", "16px", "18px"]), "10px");
              return react__WEBPACK_IMPORTED_MODULE_4___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_6__.a, {
                  sequences: sequences,
                  markerComponent: ({index: index})=>react__WEBPACK_IMPORTED_MODULE_4___default.a.createElement("div", {
                      style: {
                          width: 20,
                          display: "inline-block",
                          textAlign: "center",
                          fontSize: fontSize,
                          fontWeight: "bold",
                          color: "pink"
                      }
                  }, index)
              })
          }
          )).add("Custom styling", (function() {
              const options = {
                  sequences: sequences,
                  labelStyle: {
                      outline: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__.text)("Label style (outline)", "1px solid black")
                  },
                  markerStyle: {
                      outline: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__.text)("Marker style (outline)", "1px solid black")
                  },
                  sequenceTextColor: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__.text)("Sequence color", "blue")
              };
              return react__WEBPACK_IMPORTED_MODULE_4___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_6__.a, options)
          }
          )).add("Custom scollbars", (function() {
              const options = {
                  sequences: sequences,
                  sequenceScrollBarPositionX: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__.select)("ScrollBarPositionX", createObject(["top", "bottom"]), "top"),
                  sequenceScrollBarPositionY: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__.select)("ScrollBarPositionY", createObject(["left", "right"]), "left"),
                  sequenceOverflow: Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_8__.select)("Overflow", createObject(["scroll", "auto", "hidden"]), "scroll")
              };
              return react__WEBPACK_IMPORTED_MODULE_4___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_6__.a, options)
          }
          ))
      }
      .call(this, __webpack_require__(49)(module))
  },
  897: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.r(__webpack_exports__),
      function(module) {
          __webpack_require__(15),
          __webpack_require__(21),
          __webpack_require__(5);
          var react__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(0)
            , react__WEBPACK_IMPORTED_MODULE_3___default = __webpack_require__.n(react__WEBPACK_IMPORTED_MODULE_3__)
            , _storybook_react__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(31)
            , _lib__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(10)
            , _lib__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(318)
            , _lib__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(97)
            , _lib__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(255)
            , _lib__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(131)
            , _lib__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(254);
          function ownKeys(object, enumerableOnly) {
              var keys = Object.keys(object);
              if (Object.getOwnPropertySymbols) {
                  var symbols = Object.getOwnPropertySymbols(object);
                  enumerableOnly && (symbols = symbols.filter((function(sym) {
                      return Object.getOwnPropertyDescriptor(object, sym).enumerable
                  }
                  ))),
                  keys.push.apply(keys, symbols)
              }
              return keys
          }
          function _objectSpread(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  source = null != arguments[i] ? arguments[i] : {},
                  i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                      _defineProperty(target, key, source[key])
                  }
                  )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                      Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
                  }
                  ));
              return target
          }
          function _defineProperty(obj, key, value) {
              return key in obj ? Object.defineProperty(obj, key, {
                  value: value,
                  enumerable: !0,
                  configurable: !0,
                  writable: !0
              }) : obj[key] = value,
              obj
          }
          var withSourceLoader = __webpack_require__(38).withSource
            , __STORY__ = (__webpack_require__(38).addSource,
          '/**\n * Copyright 2018, Plotly, Inc.\n * All rights reserved.\n *\n * This source code is licensed under the MIT license found in the\n * LICENSE file in the root directory of this source tree.\n */\n\nimport React from "react";\nimport { storiesOf } from "@storybook/react";\nimport { MSAViewer } from "../lib";\n\nimport {\n  Labels,\n  OverviewBar,\n  PositionBar,\n  SequenceOverview,\n  SequenceViewer,\n} from "../lib";\n\nlet sequences = [\n  {\n    name: "seq.1",\n    sequence: "MEEPQSDPS----PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED",\n    start: 1,\n  },\n  {\n    name: "seq.2",\n    sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n    start: 40,\n  },\n  {\n    name: "seq.3",\n    sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n    start: 23,\n  },\n  {\n    name: "seq.4",\n    sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n    start: 102,\n  },\n  {\n    name: "seq.5",\n    sequence: "MEEPQSD--IEL-PLSEETFSDLWWPLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n    start: 998,\n  },\n  {\n    name: "seq.6",\n    sequence: "MEEPQEDLSSSL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n    start: 61,\n  },\n  {\n    name: "seq.7",\n    sequence: "MEEPQ---SISE-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE---LSENVAGWLEDP",\n    start: 10,\n  },\n];\n\nconst leftCoordinateStyle = {\n  textAlign: "right",\n  paddingRight: "0.25rem",\n};\n\nconst rightCoordinateStyle = {\n  paddingLeft: "0.25rem",\n};\nconst Coordinate = ({ tileHeight, style, sequence, children: coord }) => (\n  <div style={{ height: tileHeight, width: "3rem", ...style }}>\n    {sequence.start + coord}\n  </div>\n);\n\nstoriesOf("Layouting", module)\n  .add("Inverse", function () {\n    const overviewBarHeight = 50;\n    const labelsStyle = {\n      paddingTop: 20 + overviewBarHeight,\n    };\n    return <MSAViewer sequences={sequences} layout="inverse" />;\n  })\n  .add("Inverse II", function () {\n    return (\n      <MSAViewer sequences={sequences}>\n        <SequenceOverview />\n        <div style={{ display: "flex" }}>\n          <div>\n            <SequenceViewer />\n            <br />\n            <OverviewBar />\n            <PositionBar />\n          </div>\n          <Labels />\n        </div>\n      </MSAViewer>\n    );\n  })\n  .add("Full", function () {\n    return (\n      <MSAViewer sequences={sequences}>\n        <div style={{ display: "flex" }}>\n          <Labels />\n          <div>\n            <SequenceViewer />\n            <PositionBar />\n            <br />\n            <OverviewBar />\n            <br />\n            <PositionBar />\n            <OverviewBar method="information-content" />\n          </div>\n        </div>\n        <br />\n        {/* <SequenceOverview/> */}\n      </MSAViewer>\n    );\n  })\n  .add("Compact", function () {\n    return (\n      <MSAViewer sequences={sequences}>\n        <PositionBar />\n        <SequenceViewer />\n      </MSAViewer>\n    );\n  })\n  .add("Funky", function () {\n    // const options = {\n    //   sequences,\n    // };\n    // const store = createMSAStore(options);\n    const labelsStyle = {\n      paddingTop: 26,\n    };\n    return (\n      <MSAViewer sequences={sequences}>\n        <SequenceOverview />\n        <div style={{ display: "flex" }}>\n          <Labels style={labelsStyle} />\n          <div>\n            <PositionBar />\n            <SequenceViewer />\n            <PositionBar />\n            <OverviewBar />\n            <br />\n            <PositionBar />\n          </div>\n          <Labels style={labelsStyle} />\n        </div>\n        <br />\n        <SequenceOverview />\n      </MSAViewer>\n    );\n  })\n  .add("Nightingale", () => (\n    <MSAViewer sequences={sequences} layout="nightingale" height={200} />\n  ))\n  .add("Nightingale with L/R coordinates", () => (\n    <MSAViewer\n      sequences={sequences}\n      layout="nightingale"\n      leftCoordinateComponent={({ start, tileHeight, sequence }) => (\n        <Coordinate\n          tileHeight={tileHeight}\n          sequence={sequence}\n          style={leftCoordinateStyle}\n        >\n          {start}\n        </Coordinate>\n      )}\n      rightCoordinateComponent={({ end, tileHeight, sequence }) => (\n        <Coordinate\n          tileHeight={tileHeight}\n          sequence={sequence}\n          style={rightCoordinateStyle}\n        >\n          {end}\n        </Coordinate>\n      )}\n      height={200}\n    />\n  ));\n')
            , __ADDS_MAP__ = {
              "layouting--nightingale-with-l-r-coordinates": {
                  startLoc: {
                      col: 7,
                      line: 156
                  },
                  endLoc: {
                      col: 3,
                      line: 180
                  },
                  startBody: {
                      col: 43,
                      line: 156
                  },
                  endBody: {
                      col: 3,
                      line: 180
                  }
              },
              "layouting--nightingale": {
                  startLoc: {
                      col: 7,
                      line: 153
                  },
                  endLoc: {
                      col: 3,
                      line: 155
                  },
                  startBody: {
                      col: 22,
                      line: 153
                  },
                  endBody: {
                      col: 3,
                      line: 155
                  }
              },
              "layouting--funky": {
                  startLoc: {
                      col: 7,
                      line: 125
                  },
                  endLoc: {
                      col: 3,
                      line: 152
                  },
                  startBody: {
                      col: 16,
                      line: 125
                  },
                  endBody: {
                      col: 3,
                      line: 152
                  }
              },
              "layouting--compact": {
                  startLoc: {
                      col: 7,
                      line: 117
                  },
                  endLoc: {
                      col: 3,
                      line: 124
                  },
                  startBody: {
                      col: 18,
                      line: 117
                  },
                  endBody: {
                      col: 3,
                      line: 124
                  }
              },
              "layouting--full": {
                  startLoc: {
                      col: 7,
                      line: 97
                  },
                  endLoc: {
                      col: 3,
                      line: 116
                  },
                  startBody: {
                      col: 15,
                      line: 97
                  },
                  endBody: {
                      col: 3,
                      line: 116
                  }
              },
              "layouting--inverse-ii": {
                  startLoc: {
                      col: 7,
                      line: 81
                  },
                  endLoc: {
                      col: 3,
                      line: 96
                  },
                  startBody: {
                      col: 21,
                      line: 81
                  },
                  endBody: {
                      col: 3,
                      line: 96
                  }
              },
              "layouting--inverse": {
                  startLoc: {
                      col: 7,
                      line: 74
                  },
                  endLoc: {
                      col: 3,
                      line: 80
                  },
                  startBody: {
                      col: 18,
                      line: 74
                  },
                  endBody: {
                      col: 3,
                      line: 80
                  }
              }
          };
          let sequences = [{
              name: "seq.1",
              sequence: "MEEPQSDPS----PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED",
              start: 1
          }, {
              name: "seq.2",
              sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",
              start: 40
          }, {
              name: "seq.3",
              sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",
              start: 23
          }, {
              name: "seq.4",
              sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",
              start: 102
          }, {
              name: "seq.5",
              sequence: "MEEPQSD--IEL-PLSEETFSDLWWPLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",
              start: 998
          }, {
              name: "seq.6",
              sequence: "MEEPQEDLSSSL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",
              start: 61
          }, {
              name: "seq.7",
              sequence: "MEEPQ---SISE-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE---LSENVAGWLEDP",
              start: 10
          }];
          const leftCoordinateStyle = {
              textAlign: "right",
              paddingRight: "0.25rem"
          }
            , rightCoordinateStyle = {
              paddingLeft: "0.25rem"
          }
            , Coordinate = ({tileHeight: tileHeight, style: style, sequence: sequence, children: coord})=>react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("div", {
              style: _objectSpread({
                  height: tileHeight,
                  width: "3rem"
              }, style)
          }, sequence.start + coord);
          Coordinate.displayName = "Coordinate";
          var _ref = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_5__.a, {
              sequences: sequences,
              layout: "inverse"
          })
            , _ref2 = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_6__.a, null)
            , _ref3 = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("div", null, react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_7__.a, null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("br", null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_8__.a, null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_9__.a, null))
            , _ref4 = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_10__.a, null)
            , _ref5 = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_10__.a, null)
            , _ref6 = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("div", null, react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_7__.a, null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_9__.a, null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("br", null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_8__.a, null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("br", null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_9__.a, null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_8__.a, {
              method: "information-content"
          }))
            , _ref7 = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("br", null)
            , _ref8 = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_5__.a, {
              sequences: sequences
          }, react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_9__.a, null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_7__.a, null))
            , _ref9 = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_6__.a, null)
            , _ref10 = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("div", null, react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_9__.a, null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_7__.a, null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_9__.a, null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_8__.a, null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("br", null), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_9__.a, null))
            , _ref11 = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("br", null)
            , _ref12 = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_6__.a, null)
            , _ref13 = react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_5__.a, {
              sequences: sequences,
              layout: "nightingale",
              height: 200
          });
          Object(_storybook_react__WEBPACK_IMPORTED_MODULE_4__.storiesOf)("Layouting", module).addParameters({
              storySource: {
                  source: __STORY__,
                  locationsMap: __ADDS_MAP__
              }
          }).addDecorator(withSourceLoader(__STORY__, __ADDS_MAP__, "/layouts.js", [], {}, "/opt/build/repo/src/stories", {})).add("Inverse", (function() {
              return _ref
          }
          )).add("Inverse II", (function() {
              return react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_5__.a, {
                  sequences: sequences
              }, _ref2, react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("div", {
                  style: {
                      display: "flex"
                  }
              }, _ref3, _ref4))
          }
          )).add("Full", (function() {
              return react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_5__.a, {
                  sequences: sequences
              }, react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("div", {
                  style: {
                      display: "flex"
                  }
              }, _ref5, _ref6), _ref7)
          }
          )).add("Compact", (function() {
              return _ref8
          }
          )).add("Funky", (function() {
              const labelsStyle = {
                  paddingTop: 26
              };
              return react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_5__.a, {
                  sequences: sequences
              }, _ref9, react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("div", {
                  style: {
                      display: "flex"
                  }
              }, react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_10__.a, {
                  style: labelsStyle
              }), _ref10, react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_10__.a, {
                  style: labelsStyle
              })), _ref11, _ref12)
          }
          )).add("Nightingale", ()=>_ref13).add("Nightingale with L/R coordinates", ()=>react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_5__.a, {
              sequences: sequences,
              layout: "nightingale",
              leftCoordinateComponent: ({start: start, tileHeight: tileHeight, sequence: sequence})=>react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(Coordinate, {
                  tileHeight: tileHeight,
                  sequence: sequence,
                  style: leftCoordinateStyle
              }, start),
              rightCoordinateComponent: ({end: end, tileHeight: tileHeight, sequence: sequence})=>react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(Coordinate, {
                  tileHeight: tileHeight,
                  sequence: sequence,
                  style: rightCoordinateStyle
              }, end),
              height: 200
          }))
      }
      .call(this, __webpack_require__(49)(module))
  },
  898: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.r(__webpack_exports__),
      function(module) {
          __webpack_require__(15),
          __webpack_require__(30),
          __webpack_require__(24),
          __webpack_require__(7),
          __webpack_require__(21),
          __webpack_require__(5),
          __webpack_require__(26);
          var react__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(0)
            , react__WEBPACK_IMPORTED_MODULE_7___default = __webpack_require__.n(react__WEBPACK_IMPORTED_MODULE_7__)
            , _storybook_react__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(31)
            , _storybook_addon_actions__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(469)
            , _lib__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(10)
            , _lib__WEBPACK_IMPORTED_MODULE_11__ = __webpack_require__(97)
            , _lib__WEBPACK_IMPORTED_MODULE_12__ = __webpack_require__(92);
          function _extends() {
              return (_extends = Object.assign || function(target) {
                  for (var source, i = 1; i < arguments.length; i++)
                      for (var key in source = arguments[i])
                          Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
                  return target
              }
              ).apply(this, arguments)
          }
          function ownKeys(object, enumerableOnly) {
              var keys = Object.keys(object);
              if (Object.getOwnPropertySymbols) {
                  var symbols = Object.getOwnPropertySymbols(object);
                  enumerableOnly && (symbols = symbols.filter((function(sym) {
                      return Object.getOwnPropertyDescriptor(object, sym).enumerable
                  }
                  ))),
                  keys.push.apply(keys, symbols)
              }
              return keys
          }
          function _objectSpread(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  source = null != arguments[i] ? arguments[i] : {},
                  i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                      _defineProperty(target, key, source[key])
                  }
                  )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                      Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
                  }
                  ));
              return target
          }
          function _defineProperty(obj, key, value) {
              return key in obj ? Object.defineProperty(obj, key, {
                  value: value,
                  enumerable: !0,
                  configurable: !0,
                  writable: !0
              }) : obj[key] = value,
              obj
          }
          function _objectWithoutProperties(source, excluded) {
              if (null == source)
                  return {};
              var key, i, target = function _objectWithoutPropertiesLoose(source, excluded) {
                  if (null == source)
                      return {};
                  var key, i, target = {}, sourceKeys = Object.keys(source);
                  for (i = 0; i < sourceKeys.length; i++)
                      key = sourceKeys[i],
                      0 <= excluded.indexOf(key) || (target[key] = source[key]);
                  return target
              }(source, excluded);
              if (Object.getOwnPropertySymbols) {
                  var sourceSymbolKeys = Object.getOwnPropertySymbols(source);
                  for (i = 0; i < sourceSymbolKeys.length; i++)
                      key = sourceSymbolKeys[i],
                      0 <= excluded.indexOf(key) || Object.prototype.propertyIsEnumerable.call(source, key) && (target[key] = source[key])
              }
              return target
          }
          var withSourceLoader = __webpack_require__(38).withSource
            , __STORY__ = (__webpack_require__(38).addSource,
          '/**\n * Copyright 2018, Plotly, Inc.\n * All rights reserved.\n *\n * This source code is licensed under the MIT license found in the\n * LICENSE file in the root directory of this source tree.\n */\n\nimport React, { Component } from "react";\nimport { storiesOf } from "@storybook/react";\nimport { action } from "@storybook/addon-actions";\nimport { actions, MSAViewer, SequenceViewer } from "../lib";\n\nconst sequences = [\n  {\n    name: "seq.1",\n    sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED",\n  },\n  {\n    name: "seq.2",\n    sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n  {\n    name: "seq.3",\n    sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n];\n\nconst features = [\n  {\n    residues: { from: 2, to: 15 },\n    sequences: { from: 0, to: 1 },\n    id: "id-1",\n    borderColor: "blue",\n    fillColor: "black",\n  },\n  {\n    residues: { from: 25, to: 50 },\n    sequences: { from: 2, to: 2 },\n    id: "id-2",\n    borderColor: "black",\n    mouseOverBorderColor: "green",\n    fillColor: "transparent",\n    mouseOverFillColor: "transparent",\n  },\n];\n\n// storybook-action-logger doesn\'t support auto event expansion,\n// but most consoles do\nconst storyAction = (name) => {\n  const actionCallback = action(name);\n  return (e) => {\n    console.log(name, e);\n    actionCallback(e);\n  };\n};\n\nfunction Tooltip(props) {\n  const { direction, style, children, ...otherProps } = props;\n  const containerStyle = {\n    display: "inline-block",\n  };\n  const tooltipStyle = {\n    position: "relative",\n    width: "160px",\n  };\n  const textStyle = {\n    color: "#fff",\n    fontSize: "14px",\n    lineHeight: 1.2,\n    textAlign: "center",\n    backgroundColor: "#000",\n    borderRadius: "3px",\n    padding: "7px",\n  };\n  const triangleStyle = {\n    position: "absolute",\n    width: 0,\n    fontSize: 0,\n    lineHeight: 0,\n    visibility: "visible",\n    opacity: 1,\n  };\n\n  switch (direction) {\n    case "up":\n    case "down":\n      triangleStyle.borderLeft = "5px solid transparent";\n      triangleStyle.borderRight = "5px solid transparent";\n      triangleStyle.left = "50%";\n      triangleStyle.marginLeft = "-5px";\n      break;\n    case "left":\n    case "right":\n      triangleStyle.borderTop = "5px solid transparent";\n      triangleStyle.borderBottom = "5px solid transparent";\n      triangleStyle.top = "50%";\n      triangleStyle.marginTop = "-5px";\n      break;\n    default:\n  }\n\n  switch (direction) {\n    case "down":\n      triangleStyle.borderTop = "5px solid #000";\n      triangleStyle.top = "100%";\n      containerStyle.paddingBottom = "5px";\n      break;\n    case "up":\n      triangleStyle.borderBottom = "5px solid #000";\n      triangleStyle.top = "0%";\n      triangleStyle.marginTop = "-5px";\n      containerStyle.paddingTop = "5px";\n      break;\n    case "left":\n      triangleStyle.borderRight = "5px solid #000";\n      triangleStyle.marginLeft = "-5px";\n      containerStyle.paddingLeft = "5px";\n      break;\n    case "right":\n      triangleStyle.left = "100%";\n      triangleStyle.borderLeft = "5px solid #000";\n      containerStyle.paddingRight = "5px";\n      break;\n    default:\n  }\n  return (\n    <div style={{ ...containerStyle, ...style }} {...otherProps}>\n      <div style={tooltipStyle}>\n        <div style={textStyle}>{children}</div>\n        <div style={triangleStyle} />\n      </div>\n    </div>\n  );\n}\nTooltip.defaultProps = {\n  style: {},\n  direction: "down",\n};\n\nstoriesOf("Events", module)\n  .add("onResidue", () => (\n    <MSAViewer sequences={sequences}>\n      Check the console or the "Action Logger" tab for the resulting events.\n      <SequenceViewer\n        onResidueMouseEnter={storyAction("onResidueMouseEnter")}\n        onResidueMouseLeave={storyAction("onResidueMouseLeave")}\n        onResidueClick={storyAction("onResidueClick")}\n        onResidueDoubleClick={storyAction("onResidueDoubleClick")}\n      />\n    </MSAViewer>\n  ))\n  .add("onResidueClick", () => {\n    class ExtraInformation extends Component {\n      state = {};\n      onResidueClick = (e) => {\n        this.setState({ lastEvent: e });\n      };\n      render() {\n        return (\n          <div>\n            Click on a residue:\n            <MSAViewer sequences={sequences}>\n              <SequenceViewer onResidueClick={this.onResidueClick} />\n            </MSAViewer>\n            {this.state.lastEvent && (\n              <div>\n                Selection: {this.state.lastEvent.residue}\n                (from {this.state.lastEvent.sequence.name})\n              </div>\n            )}\n          </div>\n        );\n      }\n    }\n    return <ExtraInformation />;\n  })\n  .add("onFeatureClick", () => {\n    class ExtraInformation extends Component {\n      state = {};\n      onFeatureClick = (e) => {\n        this.setState({ lastEvent: e.id });\n      };\n      render() {\n        return (\n          <div>\n            Click on a feature:\n            <MSAViewer sequences={sequences}>\n              <SequenceViewer\n                onFeatureClick={this.onFeatureClick}\n                features={features}\n              />\n            </MSAViewer>\n            {this.state.lastEvent && (\n              <div>Last feature clicked: {this.state.lastEvent}</div>\n            )}\n          </div>\n        );\n      }\n    }\n    return <ExtraInformation />;\n  })\n  .add("dispatch", () => {\n    class MSADispatch extends Component {\n      onSpecificClick = (e) => {\n        this.el.updatePosition({ xPos: 100, yPos: 100 });\n      };\n      onGenericClick = (e) => {\n        const action = actions.movePosition({ xMovement: 50 });\n        this.el.dispatch(action);\n      };\n      render() {\n        return (\n          <div>\n            <MSAViewer ref={(ref) => (this.el = ref)} sequences={sequences} />\n            <button onClick={this.onSpecificClick}>Specific method</button>\n            <button onClick={this.onGenericClick}>Generic dispatch</button>\n          </div>\n        );\n      }\n    }\n    return <MSADispatch />;\n  })\n  .add("Required for Nightingale", () => {\n    class MSADispatch extends Component {\n      state = { tileWidth: 40, width: 700, aaPos: 1, highlight: null };\n      goToSpecificResidue = (aaPos) => {\n        this.setState({ aaPos });\n        // this.el.updatePositionByResidue({ aaPos });\n      };\n      modifyTileWidth = (increment) => () => {\n        // this.el.updateProp({key: \'tileWidth\', value: this.state.tileWidth + increment})\n        this.setState({ tileWidth: this.state.tileWidth + increment });\n      };\n      modifyWidth = (increment) => () => {\n        this.setState({ width: this.state.width + increment });\n      };\n      goToRegion = () => {\n        const region = {\n          start: 10,\n          end: 20,\n        };\n        this.setState({\n          tileWidth: this.state.width / (region.end + 1 - region.start),\n          aaPos: region.start,\n        });\n      };\n      highlightRegion = (n) => {\n        const highlight = {\n          sequences: {\n            from: 0,\n            to: 2,\n          },\n          residues: {\n            from: 2,\n            to: 13,\n          },\n        };\n\n        if (n === 1) this.setState({ highlight });\n        else\n          this.setState({\n            highlight: [\n              highlight,\n              {\n                ...highlight,\n                residues: {\n                  from: 20,\n                  to: 25,\n                },\n              },\n            ],\n          });\n      };\n      removeHighlightRegion = () => {\n        this.setState({ highlight: null });\n      };\n      toBeImplemented = () => console.log("Missing method");\n\n      render() {\n        const xPos = this.state.tileWidth * (this.state.aaPos - 1);\n        return (\n          <div>\n            <MSAViewer\n              ref={(ref) => (this.el = ref)}\n              sequences={sequences}\n              tileWidth={this.state.tileWidth}\n              width={this.state.width}\n              layout="compact"\n              position={{ xPos }}\n              highlight={this.state.highlight}\n            />\n            <div>\n              <button onClick={() => this.goToSpecificResidue(10.5)}>\n                GoTo Residue 10.5\n              </button>\n              <input\n                type="range"\n                min="1"\n                max="40"\n                value={this.state.aaPos}\n                onChange={(evt) => this.setState({ aaPos: evt.target.value })}\n              />\n              <code>\n                {this.state.aaPos} = {xPos}px\n              </code>\n            </div>\n            <div>\n              <b>Tile width</b>\n              <button onClick={this.modifyTileWidth(+5)}>+</button>\n              <button onClick={this.modifyTileWidth(-5)}>-</button>\n            </div>\n            <div>\n              <b>Width</b>\n              <button onClick={this.modifyWidth(+20)}>+</button>\n              <button onClick={this.modifyWidth(-20)}>-</button>\n            </div>\n            <button onClick={this.goToRegion}>\n              GoTo Residue Region [10-20]{" "}\n            </button>\n            <div>\n              <button onClick={() => this.highlightRegion(1)}>\n                Highlight Region [2-13]{" "}\n              </button>\n              <button onClick={() => this.highlightRegion(2)}>\n                Highlight Region [2-13] [20-25]{" "}\n              </button>\n              <button onClick={this.removeHighlightRegion}>\n                Remove Highlight\n              </button>\n            </div>\n            <p>Current position: [x1,x2]</p>\n          </div>\n        );\n      }\n    }\n    return <MSADispatch />;\n  })\n  .add("Tooltips (WIP)", () => {\n    class SimpleTooltip extends Component {\n      state = {};\n      onResidueMouseEnter = (e) => {\n        let direction, tooltipPosition;\n        if (e.position < 10) {\n          direction = "left";\n          tooltipPosition = {\n            top: (e.i - 0.3) * 20 + "px",\n            left: (e.position + 1) * 20 + "px",\n          };\n        } else {\n          direction = "right";\n          tooltipPosition = {\n            top: (e.i - 0.3) * 20 + "px",\n            left: e.position * 20 - 165 + "px",\n          };\n        }\n        this.setState({\n          lastEvent: e,\n          tooltipPosition,\n          direction,\n        });\n      };\n      onResidueMouseLeave = (e) => {\n        this.setState({ lastEvent: undefined });\n      };\n      render() {\n        return (\n          <div>\n            <MSAViewer sequences={sequences}>\n              <div style={{ position: "relative" }}>\n                <SequenceViewer\n                  onResidueMouseEnter={this.onResidueMouseEnter}\n                  onResidueMouseLeave={this.onResidueMouseLeave}\n                />\n                {this.state.lastEvent && (\n                  <div\n                    style={{\n                      position: "absolute",\n                      opacity: 0.8,\n                      ...this.state.tooltipPosition,\n                    }}\n                  >\n                    <Tooltip direction={this.state.direction}>\n                      {this.state.lastEvent.residue}\n                    </Tooltip>\n                  </div>\n                )}\n              </div>\n            </MSAViewer>\n          </div>\n        );\n      }\n    }\n    return <SimpleTooltip />;\n  });\n')
            , __ADDS_MAP__ = {
              "events--tooltips-wip": {
                  startLoc: {
                      col: 7,
                      line: 339
                  },
                  endLoc: {
                      col: 3,
                      line: 395
                  },
                  startBody: {
                      col: 25,
                      line: 339
                  },
                  endBody: {
                      col: 3,
                      line: 395
                  }
              },
              "events--required-for-nightingale": {
                  startLoc: {
                      col: 7,
                      line: 224
                  },
                  endLoc: {
                      col: 3,
                      line: 338
                  },
                  startBody: {
                      col: 35,
                      line: 224
                  },
                  endBody: {
                      col: 3,
                      line: 338
                  }
              },
              "events--dispatch": {
                  startLoc: {
                      col: 7,
                      line: 203
                  },
                  endLoc: {
                      col: 3,
                      line: 223
                  },
                  startBody: {
                      col: 19,
                      line: 203
                  },
                  endBody: {
                      col: 3,
                      line: 223
                  }
              },
              "events--onfeatureclick": {
                  startLoc: {
                      col: 7,
                      line: 178
                  },
                  endLoc: {
                      col: 3,
                      line: 202
                  },
                  startBody: {
                      col: 25,
                      line: 178
                  },
                  endBody: {
                      col: 3,
                      line: 202
                  }
              },
              "events--onresidueclick": {
                  startLoc: {
                      col: 7,
                      line: 153
                  },
                  endLoc: {
                      col: 3,
                      line: 177
                  },
                  startBody: {
                      col: 25,
                      line: 153
                  },
                  endBody: {
                      col: 3,
                      line: 177
                  }
              },
              "events--onresidue": {
                  startLoc: {
                      col: 7,
                      line: 142
                  },
                  endLoc: {
                      col: 3,
                      line: 152
                  },
                  startBody: {
                      col: 20,
                      line: 142
                  },
                  endBody: {
                      col: 3,
                      line: 152
                  }
              }
          };
          const sequences = [{
              name: "seq.1",
              sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED"
          }, {
              name: "seq.2",
              sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
          }, {
              name: "seq.3",
              sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
          }]
            , features = [{
              residues: {
                  from: 2,
                  to: 15
              },
              sequences: {
                  from: 0,
                  to: 1
              },
              id: "id-1",
              borderColor: "blue",
              fillColor: "black"
          }, {
              residues: {
                  from: 25,
                  to: 50
              },
              sequences: {
                  from: 2,
                  to: 2
              },
              id: "id-2",
              borderColor: "black",
              mouseOverBorderColor: "green",
              fillColor: "transparent",
              mouseOverFillColor: "transparent"
          }]
            , storyAction = name=>{
              const actionCallback = Object(_storybook_addon_actions__WEBPACK_IMPORTED_MODULE_9__.action)(name);
              return e=>{
                  console.log(name, e),
                  actionCallback(e)
              }
          }
          ;
          function Tooltip(props) {
              const direction = props.direction
                , style = props.style
                , children = props.children
                , otherProps = _objectWithoutProperties(props, ["direction", "style", "children"])
                , containerStyle = {
                  display: "inline-block"
              }
                , triangleStyle = {
                  position: "absolute",
                  width: 0,
                  fontSize: 0,
                  lineHeight: 0,
                  visibility: "visible",
                  opacity: 1
              };
              switch (direction) {
              case "up":
              case "down":
                  triangleStyle.borderLeft = "5px solid transparent",
                  triangleStyle.borderRight = "5px solid transparent",
                  triangleStyle.left = "50%",
                  triangleStyle.marginLeft = "-5px";
                  break;
              case "left":
              case "right":
                  triangleStyle.borderTop = "5px solid transparent",
                  triangleStyle.borderBottom = "5px solid transparent",
                  triangleStyle.top = "50%",
                  triangleStyle.marginTop = "-5px"
              }
              switch (direction) {
              case "down":
                  triangleStyle.borderTop = "5px solid #000",
                  triangleStyle.top = "100%",
                  containerStyle.paddingBottom = "5px";
                  break;
              case "up":
                  triangleStyle.borderBottom = "5px solid #000",
                  triangleStyle.top = "0%",
                  triangleStyle.marginTop = "-5px",
                  containerStyle.paddingTop = "5px";
                  break;
              case "left":
                  triangleStyle.borderRight = "5px solid #000",
                  triangleStyle.marginLeft = "-5px",
                  containerStyle.paddingLeft = "5px";
                  break;
              case "right":
                  triangleStyle.left = "100%",
                  triangleStyle.borderLeft = "5px solid #000",
                  containerStyle.paddingRight = "5px"
              }
              return react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", _extends({
                  style: _objectSpread({}, containerStyle, {}, style)
              }, otherProps), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", {
                  style: {
                      position: "relative",
                      width: "160px"
                  }
              }, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", {
                  style: {
                      color: "#fff",
                      fontSize: "14px",
                      lineHeight: 1.2,
                      textAlign: "center",
                      backgroundColor: "#000",
                      borderRadius: "3px",
                      padding: "7px"
                  }
              }, children), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", {
                  style: triangleStyle
              })))
          }
          Tooltip.displayName = "Tooltip",
          Tooltip.defaultProps = {
              style: {},
              direction: "down"
          };
          var _ref3 = react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("b", null, "Tile width")
            , _ref4 = react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("b", null, "Width")
            , _ref5 = react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("p", null, "Current position: [x1,x2]");
          Object(_storybook_react__WEBPACK_IMPORTED_MODULE_8__.storiesOf)("Events", module).addParameters({
              storySource: {
                  source: __STORY__,
                  locationsMap: __ADDS_MAP__
              }
          }).addDecorator(withSourceLoader(__STORY__, __ADDS_MAP__, "/events.js", [], {}, "/opt/build/repo/src/stories", {})).add("onResidue", ()=>react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_10__.a, {
              sequences: sequences
          }, 'Check the console or the "Action Logger" tab for the resulting events.', react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_11__.a, {
              onResidueMouseEnter: storyAction("onResidueMouseEnter"),
              onResidueMouseLeave: storyAction("onResidueMouseLeave"),
              onResidueClick: storyAction("onResidueClick"),
              onResidueDoubleClick: storyAction("onResidueDoubleClick")
          }))).add("onResidueClick", ()=>{
              class ExtraInformation extends react__WEBPACK_IMPORTED_MODULE_7__.Component {
                  constructor(...args) {
                      super(...args),
                      _defineProperty(this, "state", {}),
                      _defineProperty(this, "onResidueClick", e=>{
                          this.setState({
                              lastEvent: e
                          })
                      }
                      )
                  }
                  render() {
                      return react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", null, "Click on a residue:", react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_10__.a, {
                          sequences: sequences
                      }, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_11__.a, {
                          onResidueClick: this.onResidueClick
                      })), this.state.lastEvent && react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", null, "Selection: ", this.state.lastEvent.residue, "(from ", this.state.lastEvent.sequence.name, ")"))
                  }
              }
              return ExtraInformation.displayName = "ExtraInformation",
              react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(ExtraInformation, null)
          }
          ).add("onFeatureClick", ()=>{
              class ExtraInformation extends react__WEBPACK_IMPORTED_MODULE_7__.Component {
                  constructor(...args) {
                      super(...args),
                      _defineProperty(this, "state", {}),
                      _defineProperty(this, "onFeatureClick", e=>{
                          this.setState({
                              lastEvent: e.id
                          })
                      }
                      )
                  }
                  render() {
                      return react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", null, "Click on a feature:", react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_10__.a, {
                          sequences: sequences
                      }, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_11__.a, {
                          onFeatureClick: this.onFeatureClick,
                          features: features
                      })), this.state.lastEvent && react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", null, "Last feature clicked: ", this.state.lastEvent))
                  }
              }
              return ExtraInformation.displayName = "ExtraInformation",
              react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(ExtraInformation, null)
          }
          ).add("dispatch", ()=>{
              class MSADispatch extends react__WEBPACK_IMPORTED_MODULE_7__.Component {
                  constructor(...args) {
                      super(...args),
                      _defineProperty(this, "onSpecificClick", ()=>{
                          this.el.updatePosition({
                              xPos: 100,
                              yPos: 100
                          })
                      }
                      ),
                      _defineProperty(this, "onGenericClick", ()=>{
                          const action = _lib__WEBPACK_IMPORTED_MODULE_12__.a.movePosition({
                              xMovement: 50
                          });
                          this.el.dispatch(action)
                      }
                      )
                  }
                  render() {
                      return react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", null, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_10__.a, {
                          ref: _ref=>this.el = _ref,
                          sequences: sequences
                      }), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("button", {
                          onClick: this.onSpecificClick
                      }, "Specific method"), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("button", {
                          onClick: this.onGenericClick
                      }, "Generic dispatch"))
                  }
              }
              return MSADispatch.displayName = "MSADispatch",
              react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(MSADispatch, null)
          }
          ).add("Required for Nightingale", ()=>{
              class MSADispatch extends react__WEBPACK_IMPORTED_MODULE_7__.Component {
                  constructor(...args) {
                      super(...args),
                      _defineProperty(this, "state", {
                          tileWidth: 40,
                          width: 700,
                          aaPos: 1,
                          highlight: null
                      }),
                      _defineProperty(this, "goToSpecificResidue", aaPos=>{
                          this.setState({
                              aaPos: aaPos
                          })
                      }
                      ),
                      _defineProperty(this, "modifyTileWidth", increment=>()=>{
                          this.setState({
                              tileWidth: this.state.tileWidth + increment
                          })
                      }
                      ),
                      _defineProperty(this, "modifyWidth", increment=>()=>{
                          this.setState({
                              width: this.state.width + increment
                          })
                      }
                      ),
                      _defineProperty(this, "goToRegion", ()=>{
                          const region_start = 10
                            , region_end = 20;
                          this.setState({
                              tileWidth: this.state.width / (region_end + 1 - region_start),
                              aaPos: region_start
                          })
                      }
                      ),
                      _defineProperty(this, "highlightRegion", n=>{
                          const highlight = {
                              sequences: {
                                  from: 0,
                                  to: 2
                              },
                              residues: {
                                  from: 2,
                                  to: 13
                              }
                          };
                          1 === n ? this.setState({
                              highlight: highlight
                          }) : this.setState({
                              highlight: [highlight, _objectSpread({}, highlight, {
                                  residues: {
                                      from: 20,
                                      to: 25
                                  }
                              })]
                          })
                      }
                      ),
                      _defineProperty(this, "removeHighlightRegion", ()=>{
                          this.setState({
                              highlight: null
                          })
                      }
                      ),
                      _defineProperty(this, "toBeImplemented", ()=>console.log("Missing method"))
                  }
                  render() {
                      const xPos = this.state.tileWidth * (this.state.aaPos - 1);
                      return react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", null, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_10__.a, {
                          ref: _ref2=>this.el = _ref2,
                          sequences: sequences,
                          tileWidth: this.state.tileWidth,
                          width: this.state.width,
                          layout: "compact",
                          position: {
                              xPos: xPos
                          },
                          highlight: this.state.highlight
                      }), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", null, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("button", {
                          onClick: ()=>this.goToSpecificResidue(10.5)
                      }, "GoTo Residue 10.5"), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("input", {
                          type: "range",
                          min: "1",
                          max: "40",
                          value: this.state.aaPos,
                          onChange: evt=>this.setState({
                              aaPos: evt.target.value
                          })
                      }), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("code", null, this.state.aaPos, " = ", xPos, "px")), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", null, _ref3, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("button", {
                          onClick: this.modifyTileWidth(5)
                      }, "+"), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("button", {
                          onClick: this.modifyTileWidth(-5)
                      }, "-")), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", null, _ref4, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("button", {
                          onClick: this.modifyWidth(20)
                      }, "+"), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("button", {
                          onClick: this.modifyWidth(-20)
                      }, "-")), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("button", {
                          onClick: this.goToRegion
                      }, "GoTo Residue Region [10-20]", " "), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", null, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("button", {
                          onClick: ()=>this.highlightRegion(1)
                      }, "Highlight Region [2-13]", " "), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("button", {
                          onClick: ()=>this.highlightRegion(2)
                      }, "Highlight Region [2-13] [20-25]", " "), react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("button", {
                          onClick: this.removeHighlightRegion
                      }, "Remove Highlight")), _ref5)
                  }
              }
              return MSADispatch.displayName = "MSADispatch",
              react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(MSADispatch, null)
          }
          ).add("Tooltips (WIP)", ()=>{
              class SimpleTooltip extends react__WEBPACK_IMPORTED_MODULE_7__.Component {
                  constructor(...args) {
                      super(...args),
                      _defineProperty(this, "state", {}),
                      _defineProperty(this, "onResidueMouseEnter", e=>{
                          let direction, tooltipPosition;
                          10 > e.position ? (direction = "left",
                          tooltipPosition = {
                              top: 20 * (e.i - .3) + "px",
                              left: 20 * (e.position + 1) + "px"
                          }) : (direction = "right",
                          tooltipPosition = {
                              top: 20 * (e.i - .3) + "px",
                              left: 20 * e.position - 165 + "px"
                          }),
                          this.setState({
                              lastEvent: e,
                              tooltipPosition: tooltipPosition,
                              direction: direction
                          })
                      }
                      ),
                      _defineProperty(this, "onResidueMouseLeave", ()=>{
                          this.setState({
                              lastEvent: void 0
                          })
                      }
                      )
                  }
                  render() {
                      return react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", null, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_10__.a, {
                          sequences: sequences
                      }, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", {
                          style: {
                              position: "relative"
                          }
                      }, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_11__.a, {
                          onResidueMouseEnter: this.onResidueMouseEnter,
                          onResidueMouseLeave: this.onResidueMouseLeave
                      }), this.state.lastEvent && react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement("div", {
                          style: _objectSpread({
                              position: "absolute",
                              opacity: .8
                          }, this.state.tooltipPosition)
                      }, react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(Tooltip, {
                          direction: this.state.direction
                      }, this.state.lastEvent.residue)))))
                  }
              }
              return SimpleTooltip.displayName = "SimpleTooltip",
              react__WEBPACK_IMPORTED_MODULE_7___default.a.createElement(SimpleTooltip, null)
          }
          )
      }
      .call(this, __webpack_require__(49)(module))
  },
  911: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.r(__webpack_exports__),
      function(module) {
          __webpack_require__(24),
          __webpack_require__(7),
          __webpack_require__(26);
          var react__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(0)
            , react__WEBPACK_IMPORTED_MODULE_3___default = __webpack_require__.n(react__WEBPACK_IMPORTED_MODULE_3__)
            , _storybook_react__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(31)
            , _lib__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(10)
            , _storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(9)
            , lodash_es__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(317)
            , lodash_es__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(917);
          function _extends() {
              return (_extends = Object.assign || function(target) {
                  for (var source, i = 1; i < arguments.length; i++)
                      for (var key in source = arguments[i])
                          Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
                  return target
              }
              ).apply(this, arguments)
          }
          function _defineProperty(obj, key, value) {
              return key in obj ? Object.defineProperty(obj, key, {
                  value: value,
                  enumerable: !0,
                  configurable: !0,
                  writable: !0
              }) : obj[key] = value,
              obj
          }
          var withSourceLoader = __webpack_require__(38).withSource
            , __STORY__ = (__webpack_require__(38).addSource,
          '/**\n * Copyright 2018, Plotly, Inc.\n * All rights reserved.\n *\n * This source code is licensed under the MIT license found in the\n * LICENSE file in the root directory of this source tree.\n */\n\nimport React, { Component } from "react";\nimport { storiesOf } from "@storybook/react";\nimport { MSAViewer } from "../lib";\nimport { number, button, withKnobs } from "@storybook/addon-knobs";\nimport { repeat, times } from "lodash-es";\n\nconst rawSequences = [\n  "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED",\n  "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n];\n\nconst sequences = [];\n\ntimes(1000, (i) => {\n  sequences.push({\n    name: `sequence ${i}`,\n    sequence: repeat(rawSequences[i % 3], 5),\n  });\n});\n\nstoriesOf("Navigation", module)\n  .addDecorator(withKnobs)\n  .add("Move viewpoint", function () {\n    const options = {\n      sequences,\n      sequenceOverflow: "hidden",\n      width: 600,\n      height: 600,\n    };\n    const mov = number("Movement by", 10, {\n      range: true,\n      min: 1,\n      max: 100,\n      step: 1,\n    });\n\n    class MSA extends Component {\n      state = {\n        position: {\n          xPos: 100,\n          yPos: 100,\n        },\n      };\n      moveLeft = () => this.move({ x: -mov });\n      moveRight = () => this.move({ x: mov });\n      moveTop = () => this.move({ y: -mov });\n      moveBottom = () => this.move({ y: mov });\n\n      move = ({ x = 0, y = 0 }) => {\n        const { xPos, yPos } = this.state.position;\n        const position = {\n          xPos: xPos + x,\n          yPos: yPos + y,\n        };\n        this.setState({ position });\n      };\n      startLoop = () => {\n        if (this.frame) return;\n        this.counter = 0;\n        this.loop = () => {\n          const xMov = this.counter % 40 < 20 ? 5 : -5;\n          this.counter++;\n          this.move({ x: xMov, y: xMov });\n          this.frame = window.requestAnimationFrame(this.loop);\n        };\n        this.loop();\n      };\n      endLoop = () => {\n        window.cancelAnimationFrame(this.frame);\n      };\n      render() {\n        return (\n          <div>\n            <MSAViewer {...this.props} {...this.state} />\n            <div>\n              <button onClick={this.moveLeft}>Left</button>\n              <button onClick={this.moveRight}>Right</button>\n              <button onClick={this.moveTop}>Top</button>\n              <button onClick={this.moveBottom}>Bottom</button>\n            </div>\n            <button onClick={this.startLoop}>Start Loop</button>\n            <button onClick={this.endLoop}>End Loop</button>\n          </div>\n        );\n      }\n    }\n    return <MSA {...options} />;\n  });\n')
            , __ADDS_MAP__ = {
              "navigation--move-viewpoint": {
                  startLoc: {
                      col: 7,
                      line: 32
                  },
                  endLoc: {
                      col: 3,
                      line: 97
                  },
                  startBody: {
                      col: 25,
                      line: 32
                  },
                  endBody: {
                      col: 3,
                      line: 97
                  }
              }
          };
          const rawSequences = ["MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED", "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP", "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"]
            , sequences = [];
          Object(lodash_es__WEBPACK_IMPORTED_MODULE_7__.a)(1e3, i=>{
              sequences.push({
                  name: "sequence ".concat(i),
                  sequence: Object(lodash_es__WEBPACK_IMPORTED_MODULE_8__.a)(rawSequences[i % 3], 5)
              })
          }
          ),
          Object(_storybook_react__WEBPACK_IMPORTED_MODULE_4__.storiesOf)("Navigation", module).addParameters({
              storySource: {
                  source: __STORY__,
                  locationsMap: __ADDS_MAP__
              }
          }).addDecorator(withSourceLoader(__STORY__, __ADDS_MAP__, "/navigation.js", [], {}, "/opt/build/repo/src/stories", {})).addDecorator(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_6__.withKnobs).add("Move viewpoint", (function() {
              const mov = Object(_storybook_addon_knobs__WEBPACK_IMPORTED_MODULE_6__.number)("Movement by", 10, {
                  range: !0,
                  min: 1,
                  max: 100,
                  step: 1
              });
              class MSA extends react__WEBPACK_IMPORTED_MODULE_3__.Component {
                  constructor(...args) {
                      super(...args),
                      _defineProperty(this, "state", {
                          position: {
                              xPos: 100,
                              yPos: 100
                          }
                      }),
                      _defineProperty(this, "moveLeft", ()=>this.move({
                          x: -mov
                      })),
                      _defineProperty(this, "moveRight", ()=>this.move({
                          x: mov
                      })),
                      _defineProperty(this, "moveTop", ()=>this.move({
                          y: -mov
                      })),
                      _defineProperty(this, "moveBottom", ()=>this.move({
                          y: mov
                      })),
                      _defineProperty(this, "move", ({x: x=0, y: y=0})=>{
                          const _this$state$position = this.state.position
                            , xPos = _this$state$position.xPos
                            , yPos = _this$state$position.yPos;
                          this.setState({
                              position: {
                                  xPos: xPos + x,
                                  yPos: yPos + y
                              }
                          })
                      }
                      ),
                      _defineProperty(this, "startLoop", ()=>{
                          this.frame || (this.counter = 0,
                          this.loop = ()=>{
                              const xMov = 20 > this.counter % 40 ? 5 : -5;
                              this.counter++,
                              this.move({
                                  x: xMov,
                                  y: xMov
                              }),
                              this.frame = window.requestAnimationFrame(this.loop)
                          }
                          ,
                          this.loop())
                      }
                      ),
                      _defineProperty(this, "endLoop", ()=>{
                          window.cancelAnimationFrame(this.frame)
                      }
                      )
                  }
                  render() {
                      return react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("div", null, react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_5__.a, _extends({}, this.props, this.state)), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("div", null, react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("button", {
                          onClick: this.moveLeft
                      }, "Left"), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("button", {
                          onClick: this.moveRight
                      }, "Right"), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("button", {
                          onClick: this.moveTop
                      }, "Top"), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("button", {
                          onClick: this.moveBottom
                      }, "Bottom")), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("button", {
                          onClick: this.startLoop
                      }, "Start Loop"), react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement("button", {
                          onClick: this.endLoop
                      }, "End Loop"))
                  }
              }
              return MSA.displayName = "MSA",
              react__WEBPACK_IMPORTED_MODULE_3___default.a.createElement(MSA, {
                  sequences: sequences,
                  sequenceOverflow: "hidden",
                  width: 600,
                  height: 600
              })
          }
          ))
      }
      .call(this, __webpack_require__(49)(module))
  },
  912: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.r(__webpack_exports__),
      function(module) {
          var react__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(0)
            , react__WEBPACK_IMPORTED_MODULE_0___default = __webpack_require__.n(react__WEBPACK_IMPORTED_MODULE_0__)
            , _storybook_react__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(31)
            , _lib__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(97)
            , _lib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(37)
            , _lib__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(16)
            , _lib__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(10)
            , withSourceLoader = __webpack_require__(38).withSource
            , __STORY__ = (__webpack_require__(38).addSource,
          '/**\n * Copyright 2018, Plotly, Inc.\n * All rights reserved.\n *\n * This source code is licensed under the MIT license found in the\n * LICENSE file in the root directory of this source tree.\n */\n\nimport React, { Component } from "react";\nimport { storiesOf } from "@storybook/react";\nimport {\n  msaConnect,\n  MSAViewer,\n  withPositionStore,\n  SequenceViewer,\n} from "../lib";\n\nconst sequences = [\n  {\n    name: "seq.1",\n    sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED",\n  },\n  {\n    name: "seq.2",\n    sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n  {\n    name: "seq.3",\n    sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP",\n  },\n];\n\nstoriesOf("Plugins", module).add("My first plugin", function () {\n  class MyFirstMSAPluginComponent extends Component {\n    // called on every position update (e.g. mouse movement or scrolling)\n    shouldRerender(newPosition) {\n      return true;\n    }\n    render() {\n      return (\n        <div>\n          x: {this.props.position.xPos}, y: {this.props.position.yPos}\n        </div>\n      );\n    }\n  }\n\n  // inject position awareness (this is done to avoid react tree computations)\n  // "performance is the root of all evil"\n  const MyFirstMSAPluginConnected = withPositionStore(\n    MyFirstMSAPluginComponent\n  );\n\n  // select attributes from the main redux store\n  const mapStateToProps = (state) => {\n    return {\n      height: state.props.height,\n      sequences: state.sequences,\n    };\n  };\n\n  // subscribe to the main redux store\n  const MyFirstMSAPlugin = msaConnect(mapStateToProps)(\n    MyFirstMSAPluginConnected\n  );\n\n  return (\n    <MSAViewer sequences={sequences} height={60}>\n      <SequenceViewer />\n      <MyFirstMSAPlugin />\n    </MSAViewer>\n  );\n});\n')
            , __ADDS_MAP__ = {
              "plugins--my-first-plugin": {
                  startLoc: {
                      col: 33,
                      line: 33
                  },
                  endLoc: {
                      col: 1,
                      line: 73
                  },
                  startBody: {
                      col: 52,
                      line: 33
                  },
                  endBody: {
                      col: 1,
                      line: 73
                  }
              }
          };
          const sequences = [{
              name: "seq.1",
              sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED"
          }, {
              name: "seq.2",
              sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
          }, {
              name: "seq.3",
              sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
          }];
          var _ref = react__WEBPACK_IMPORTED_MODULE_0___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_2__.a, null);
          Object(_storybook_react__WEBPACK_IMPORTED_MODULE_1__.storiesOf)("Plugins", module).addParameters({
              storySource: {
                  source: __STORY__,
                  locationsMap: __ADDS_MAP__
              }
          }).addDecorator(withSourceLoader(__STORY__, __ADDS_MAP__, "/plugins.js", [], {}, "/opt/build/repo/src/stories", {})).add("My first plugin", (function() {
              class MyFirstMSAPluginComponent extends react__WEBPACK_IMPORTED_MODULE_0__.Component {
                  shouldRerender() {
                      return !0
                  }
                  render() {
                      return react__WEBPACK_IMPORTED_MODULE_0___default.a.createElement("div", null, "x: ", this.props.position.xPos, ", y: ", this.props.position.yPos)
                  }
              }
              MyFirstMSAPluginComponent.displayName = "MyFirstMSAPluginComponent";
              const MyFirstMSAPluginConnected = Object(_lib__WEBPACK_IMPORTED_MODULE_3__.a)(MyFirstMSAPluginComponent)
                , MyFirstMSAPlugin = Object(_lib__WEBPACK_IMPORTED_MODULE_4__.a)(state=>({
                  height: state.props.height,
                  sequences: state.sequences
              }))(MyFirstMSAPluginConnected);
              return react__WEBPACK_IMPORTED_MODULE_0___default.a.createElement(_lib__WEBPACK_IMPORTED_MODULE_5__.a, {
                  sequences: sequences,
                  height: 60
              }, _ref, react__WEBPACK_IMPORTED_MODULE_0___default.a.createElement(MyFirstMSAPlugin, null))
          }
          ))
      }
      .call(this, __webpack_require__(49)(module))
  },
  92: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__.d(__webpack_exports__, "a", (function() {
          return actions
      }
      ));
      __webpack_require__(15),
      __webpack_require__(21),
      __webpack_require__(5),
      __webpack_require__(127);
      var _store_positionReducers__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(47)
        , _components_MSAViewer__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(10);
      function ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function _defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      const actions = function _objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                  _defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }({}, __webpack_require__(8).b, {}, _store_positionReducers__WEBPACK_IMPORTED_MODULE_4__.a);
      _components_MSAViewer__WEBPACK_IMPORTED_MODULE_5__.a
  },
  93: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      function shallowEqual(objA, objB, compare, compareContext) {
          var ret = compare ? compare.call(compareContext, objA, objB) : void 0;
          if (void 0 !== ret)
              return !!ret;
          if (objA === objB)
              return !0;
          if ("object" != typeof objA || !objA || "object" != typeof objB || !objB)
              return !1;
          var keysA = Object.keys(objA)
            , keysB = Object.keys(objB);
          if (keysA.length !== keysB.length)
              return !1;
          for (var key, bHasOwnProperty = Object.prototype.hasOwnProperty.bind(objB), idx = 0; idx < keysA.length; idx++) {
              if (!bHasOwnProperty(key = keysA[idx]))
                  return !1;
              var valueA = objA[key]
                , valueB = objB[key];
              if (!1 === (ret = compare ? compare.call(compareContext, valueA, valueB, key) : void 0) || void 0 === ret && valueA !== valueB)
                  return !1
          }
          return !0
      }
      __webpack_require__.d(__webpack_exports__, "a", (function() {
          return shallowEqual
      }
      ))
  },
  97: function(module, __webpack_exports__, __webpack_require__) {
      "use strict";
      __webpack_require__(45),
      __webpack_require__(33),
      __webpack_require__(15),
      __webpack_require__(103),
      __webpack_require__(24),
      __webpack_require__(35),
      __webpack_require__(178),
      __webpack_require__(21),
      __webpack_require__(52),
      __webpack_require__(5),
      __webpack_require__(26);
      var prop_types = __webpack_require__(1)
        , prop_types_default = __webpack_require__.n(prop_types)
        , pick = __webpack_require__(319)
        , isEqual = __webpack_require__(922)
        , clamp = __webpack_require__(188)
        , floor = __webpack_require__(320)
        , react = (__webpack_require__(7),
      __webpack_require__(0))
        , react_default = __webpack_require__.n(react)
        , omit = __webpack_require__(924);
      function _defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      class Mouse {
      }
      _defineProperty(Mouse, "rel", (function(e) {
          if (void 0 !== e.changedTouches)
              return Mouse.rel(e.changedTouches[0]);
          let mouseX = e.offsetX
            , mouseY = e.offsetY;
          if (void 0 === mouseX) {
              const target = e.target || e.srcElement
                , rect = target.getBoundingClientRect();
              if (mouseX = e.clientX - rect.left,
              mouseY = e.clientY - rect.top,
              void 0 === mouseX && (mouseX = e.pageX - target.offsetLeft,
              mouseY = e.pageY - target.offsetTop,
              void 0 === mouseX))
                  return void console.log(e, "No mouse event defined.")
          }
          return [mouseX, mouseY]
      }
      )),
      _defineProperty(Mouse, "abs", (function(e) {
          if (void 0 !== e.changedTouches)
              return Mouse.abs(e.changedTouches[0]);
          let mouseX = e.pageX
            , mouseY = e.pageY;
          return void 0 === mouseX && (mouseX = e.layerX,
          mouseY = e.layerY),
          void 0 === mouseX && (mouseX = e.clientX,
          mouseY = e.clientY),
          void 0 === mouseX && (mouseX = e.x,
          mouseY = e.y),
          [mouseX, mouseY]
      }
      )),
      _defineProperty(Mouse, "wheelDelta", (function(e) {
          let delta = [e.deltaX, e.deltaY];
          return void 0 === delta[0] && e.mozMovementX && (delta = [0, e.mozMovementX]),
          isNaN(delta[0]) && (delta[0] = 0),
          isNaN(delta[1]) && (delta[1] = 0),
          delta
      }
      ));
      var mouse = Mouse
        , AutoscaleIcon = __webpack_require__(466)
        , AutoscaleIcon_default = __webpack_require__.n(AutoscaleIcon)
        , ZoomPlusIcon = __webpack_require__(467)
        , ZoomPlusIcon_default = __webpack_require__.n(ZoomPlusIcon)
        , ZoomMinusIcon = __webpack_require__(468)
        , ZoomMinusIcon_default = __webpack_require__.n(ZoomMinusIcon);
      function ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function ModBar_defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      var _ref = react_default.a.createElement("svg", {
          height: "1em",
          width: "1.542em",
          viewBox: "0 0 1542 1000"
      }, react_default.a.createElement("path", {
          d: "m0-10h182v-140h-182v140z m228 146h183v-286h-183v286z m225 714h182v-1000h-182v1000z m225-285h182v-715h-182v715z m225 142h183v-857h-183v857z m231-428h182v-429h-182v429z m225-291h183v-138h-183v138z",
          transform: "matrix(1 0 0 -1 0 850)",
          fill: "rgb(68, 122, 219)"
      }));
      function PlotlyIcon() {
          return _ref
      }
      PlotlyIcon.displayName = "PlotlyIcon";
      class ModBar_ModBar extends react.PureComponent {
          render() {
              const style = function _objectSpread(target) {
                  for (var source, i = 1; i < arguments.length; i++)
                      source = null != arguments[i] ? arguments[i] : {},
                      i % 2 ? ownKeys(Object(source), !0).forEach((function(key) {
                          ModBar_defineProperty(target, key, source[key])
                      }
                      )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : ownKeys(Object(source)).forEach((function(key) {
                          Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
                      }
                      ));
                  return target
              }({
                  opacity: .9,
                  backgroundColor: "white"
              }, this.props.style)
                , linkStyle = {
                  position: "relative",
                  fontSize: "16px",
                  padding: "3px 4px",
                  cursor: "pointer",
                  lineHeight: "normal",
                  textDecoration: "none",
                  color: "black"
              };
              return react_default.a.createElement("div", {
                  style: style
              }, react_default.a.createElement("div", {
                  style: linkStyle
              }, react_default.a.createElement(ZoomPlusIcon_default.a, {
                  width: 20,
                  height: 20
              })), react_default.a.createElement("div", {
                  style: linkStyle
              }, react_default.a.createElement(ZoomMinusIcon_default.a, {
                  width: 20,
                  height: 20
              })), react_default.a.createElement("div", {
                  style: linkStyle
              }, react_default.a.createElement(AutoscaleIcon_default.a, {
                  width: 20,
                  height: 20
              })), react_default.a.createElement("a", {
                  href: "https://plot.ly/",
                  target: "_blank",
                  rel: "noreferrer noopener",
                  "data-title": "Produced with Plotly",
                  style: linkStyle
              }, react_default.a.createElement(PlotlyIcon, {
                  width: 20,
                  height: 20
              })))
          }
      }
      ModBar_ModBar.displayName = "ModBar",
      ModBar_ModBar.__docgenInfo = {
          description: "",
          methods: [],
          displayName: "ModBar"
      };
      var components_ModBar = ModBar_ModBar;
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/ModBar.js"] = {
          name: "ModBar",
          docgenInfo: ModBar_ModBar.__docgenInfo,
          path: "src/components/ModBar.js"
      });
      var positionReducers = __webpack_require__(47)
        , withPositionStore = __webpack_require__(37)
        , requestAnimation = __webpack_require__(153);
      function FakeScroll_defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      var FakeScroll_ref = react_default.a.createElement("div", null);
      class FakeScroll_FakeScroll extends react.PureComponent {
          constructor(props) {
              super(props),
              FakeScroll_defineProperty(this, "onScroll", ()=>{
                  Object(requestAnimation.a)(this, ()=>{
                      const _this$shouldShow = this.shouldShow()
                        , showX = _this$shouldShow.showX
                        , showY = _this$shouldShow.showY
                        , movement = {
                          xMovement: this.el.current.scrollLeft - (showX ? this.props.position.xPos : 0),
                          yMovement: this.el.current.scrollTop - (showY ? this.props.position.yPos : 0)
                      };
                      this.props.positionDispatch(Object(positionReducers.b)(movement))
                  }
                  )
              }
              ),
              FakeScroll_defineProperty(this, "updateScrollPosition", ()=>{
                  this.el && this.el.current && (this.el.current.scrollTop = this.props.position.yPos,
                  this.el.current.scrollLeft = this.props.position.xPos)
              }
              ),
              this.el = react_default.a.createRef()
          }
          checkOverflow(overflow, {withX: withX=!1, withY: withY=!1}) {
              let show = !1;
              switch (overflow) {
              case "auto":
                  withX && (show |= this.props.fullWidth > this.props.width),
                  withY && (show |= this.props.fullHeight > this.props.height);
                  break;
              case "hidden":
                  show = !1;
                  break;
              case "scroll":
                  show = !0
              }
              return show
          }
          shouldShow() {
              const withX = {
                  withX: !0
              }
                , withY = {
                  withY: !0
              }
                , overflowX = "initial" === this.props.overflowX ? this.props.overflow : this.props.overflowX
                , overflowY = "initial" === this.props.overflowY ? this.props.overflow : this.props.overflowY;
              return {
                  showX: this.checkOverflow(overflowX, withX) && this.checkOverflow(this.props.overflow, withX),
                  showY: this.checkOverflow(overflowY, withY) && this.checkOverflow(this.props.overflow, withY)
              }
          }
          render() {
              const _this$props = this.props
                , width = _this$props.width
                , height = _this$props.height
                , fullWidth = _this$props.fullWidth
                , fullHeight = _this$props.fullHeight
                , style = {
                  position: "absolute",
                  overflowX: "auto",
                  overflowY: "auto",
                  width: width,
                  height: height,
                  transform: ""
              }
                , _this$shouldShow2 = this.shouldShow()
                , showX = _this$shouldShow2.showX
                , showY = _this$shouldShow2.showY
                , childStyle = {
                  height: 1,
                  width: 1
              };
              return showY || showX ? (showX && (childStyle.width = fullWidth,
              style.overflowX = "scroll",
              "top" === this.props.positionX && (style.transform += "rotateX(180deg)")),
              showY && (childStyle.height = fullHeight,
              style.overflowY = "scroll",
              "left" === this.props.positionY && (style.transform += "rotateY(180deg)")),
              react_default.a.createElement("div", {
                  style: style,
                  onScroll: this.onScroll,
                  ref: this.el
              }, react_default.a.createElement("div", {
                  style: childStyle
              }))) : FakeScroll_ref
          }
      }
      FakeScroll_FakeScroll.displayName = "FakeScroll",
      FakeScroll_FakeScroll.defaultProps = {
          overflow: "auto",
          overflowX: "initial",
          overflowY: "initial",
          positionX: "bottom",
          positionY: "right",
          scrollBarWidth: 5
      },
      FakeScroll_FakeScroll.propTypes = {
          overflow: prop_types_default.a.oneOf(["hidden", "auto", "scroll"]),
          overflowX: prop_types_default.a.oneOf(["hidden", "auto", "scroll", "initial"]),
          overflowY: prop_types_default.a.oneOf(["hidden", "auto", "scroll", "initial"]),
          width: prop_types_default.a.number.isRequired,
          height: prop_types_default.a.number.isRequired,
          fullHeight: prop_types_default.a.number.isRequired,
          fullWidth: prop_types_default.a.number.isRequired,
          positionX: prop_types_default.a.oneOf(["top", "bottom"]),
          positionY: prop_types_default.a.oneOf(["left", "right"])
      },
      FakeScroll_FakeScroll.__docgenInfo = {
          description: "Creates a DOM element with absolute position that can have scrollbars.\nHowever, no actual content is displayed by this element.",
          methods: [{
              name: "onScroll",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "e",
                  type: null
              }],
              returns: null
          }, {
              name: "updateScrollPosition",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "checkOverflow",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "overflow",
                  type: null
              }, {
                  name: "{ withX = false, withY = false }",
                  type: null
              }],
              returns: null
          }, {
              name: "shouldShow",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }],
          displayName: "FakeScroll",
          props: {
              overflow: {
                  defaultValue: {
                      value: '"auto"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"hidden"',
                          computed: !1
                      }, {
                          value: '"auto"',
                          computed: !1
                      }, {
                          value: '"scroll"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: ""
              },
              overflowX: {
                  defaultValue: {
                      value: '"initial"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"hidden"',
                          computed: !1
                      }, {
                          value: '"auto"',
                          computed: !1
                      }, {
                          value: '"scroll"',
                          computed: !1
                      }, {
                          value: '"initial"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: ""
              },
              overflowY: {
                  defaultValue: {
                      value: '"initial"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"hidden"',
                          computed: !1
                      }, {
                          value: '"auto"',
                          computed: !1
                      }, {
                          value: '"scroll"',
                          computed: !1
                      }, {
                          value: '"initial"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: ""
              },
              positionX: {
                  defaultValue: {
                      value: '"bottom"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"top"',
                          computed: !1
                      }, {
                          value: '"bottom"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: ""
              },
              positionY: {
                  defaultValue: {
                      value: '"right"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"left"',
                          computed: !1
                      }, {
                          value: '"right"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: ""
              },
              scrollBarWidth: {
                  defaultValue: {
                      value: "5",
                      computed: !1
                  },
                  required: !1
              },
              width: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              height: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              fullHeight: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              fullWidth: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              }
          }
      };
      var Canvas_FakeScroll = Object(withPositionStore.a)(FakeScroll_FakeScroll, {
          withX: !0,
          withY: !0
      });
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/Canvas/FakeScroll.js"] = {
          name: "FakeScroll",
          docgenInfo: FakeScroll_FakeScroll.__docgenInfo,
          path: "src/components/Canvas/FakeScroll.js"
      });
      var autobind = __webpack_require__(77);
      function _extends() {
          return (_extends = Object.assign || function(target) {
              for (var source, i = 1; i < arguments.length; i++)
                  for (var key in source = arguments[i])
                      Object.prototype.hasOwnProperty.call(source, key) && (target[key] = source[key]);
              return target
          }
          ).apply(this, arguments)
      }
      function DraggingComponent_ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function DraggingComponent_objectSpread(target) {
          for (var source, i = 1; i < arguments.length; i++)
              source = null != arguments[i] ? arguments[i] : {},
              i % 2 ? DraggingComponent_ownKeys(Object(source), !0).forEach((function(key) {
                  DraggingComponent_defineProperty(target, key, source[key])
              }
              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : DraggingComponent_ownKeys(Object(source)).forEach((function(key) {
                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
              }
              ));
          return target
      }
      function DraggingComponent_defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      class DraggingComponent_DraggingComponent extends react.PureComponent {
          constructor(props) {
              super(props),
              DraggingComponent_defineProperty(this, "state", {
                  mouse: {
                      isMouseWithin: !1,
                      cursorState: "grab"
                  }
              }),
              DraggingComponent_defineProperty(this, "currentContext", 1),
              this.canvasBuffers = [react_default.a.createRef(), react_default.a.createRef()],
              this.container = react_default.a.createRef(),
              Object(autobind.a)(this, "onMouseEnter", "onMouseLeave", "onMouseDown", "onMouseUp", "onMouseMove", "onTouchStart", "onTouchMove", "onTouchEnd", "onTouchCancel", "onClick", "onDoubleClick", "draw"),
              this.onViewpointChange(),
              this.mouseMovePosition = void 0,
              this.touchMovePosition = void 0,
              this.isInDragPhase = void 0,
              this.mouseHasMoved = void 0
          }
          drawScene() {
              console.warn("drawScene is unimplemented.")
          }
          onPositionUpdate() {
              console.warn("onPositionUpdate is unimplemented.")
          }
          onViewpointChange() {}
          componentDidMount() {
              this.ctxBuffers = [this.canvasBuffers[0].current.getContext("2d", {
                  alpha: "false"
              }), this.canvasBuffers[1].current.getContext("2d", {
                  alpha: "false"
              })],
              this.swapContexts(),
              this.container.current.addEventListener("mouseenter", this.onMouseEnter),
              this.container.current.addEventListener("mouseleave", this.onMouseLeave),
              this.container.current.addEventListener("touchstart", this.onTouchStart),
              this.container.current.addEventListener("touchmove", this.onTouchMove),
              this.container.current.addEventListener("touchend", this.onTouchEnd),
              this.container.current.addEventListener("touchcancel", this.onTouchCancel),
              this.container.current.addEventListener("click", this.onClick),
              this.container.current.addEventListener("dblclick", this.onDoubleClick),
              this.props.sequenceDisableDragging || (this.container.current.addEventListener("mousedown", this.onMouseDown),
              this.container.current.addEventListener("mouseup", this.onMouseUp),
              this.container.current.addEventListener("mousemove", this.onMouseMove))
          }
          swapContexts() {
              const current = this.currentContext;
              this.canvasBuffers[current].current.style.visibility = "visible";
              const next = (this.currentContext + 1) % 2;
              this.canvasBuffers[next].current.style.visibility = "hidden",
              this.currentContext = next,
              this.ctx = this.ctxBuffers[next]
          }
          draw() {
              this.ctx && (this.ctx.clearRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height),
              this.onViewpointChange(),
              this.drawScene(),
              this.swapContexts())
          }
          onClick() {}
          onDoubleClick() {}
          onMouseDown(e) {
              this.startDragPhase(e)
          }
          onMouseMove(e) {
              if (void 0 === this.isInDragPhase)
                  return;
              this.mouseHasMoved = !0;
              const pos = mouse.abs(e)
                , oldPos = this.mouseMovePosition;
              Object(requestAnimation.a)(this, ()=>{
                  this.mouseMovePosition = pos,
                  this.onPositionUpdate(oldPos, this.mouseMovePosition)
              }
              )
          }
          onMouseUp() {
              this.stopDragPhase()
          }
          onMouseEnter() {
              this.setState(prevState=>({
                  mouse: DraggingComponent_objectSpread({}, prevState.mouse, {
                      isMouseWithin: !0
                  })
              }))
          }
          onMouseLeave() {
              this.stopHoverPhase(),
              this.stopDragPhase()
          }
          onTouchStart(e) {
              console.log("touchstart", e),
              this.startDragPhase(e)
          }
          onTouchMove(e) {
              void 0 !== this.isInDragPhase && (console.log("touchmove", e),
              this.onMouseMove(e))
          }
          onTouchEnd() {
              this.stopDragPhase()
          }
          onTouchCancel() {
              this.stopDragPhase()
          }
          startDragPhase(e) {
              this.mouseMovePosition = mouse.abs(e),
              this.mouseHasMoved = void 0,
              this.isInDragPhase = !0,
              this.setState(prevState=>({
                  mouse: DraggingComponent_objectSpread({}, prevState.mouse, {
                      cursorState: "grabbing"
                  })
              }))
          }
          stopHoverPhase() {
              this.setState(prevState=>({
                  mouse: DraggingComponent_objectSpread({}, prevState.mouse, {
                      isMouseWithin: !1
                  })
              }))
          }
          stopDragPhase() {
              this.isInDragPhase = void 0,
              this.setState(prevState=>({
                  mouse: DraggingComponent_objectSpread({}, prevState.mouse, {
                      cursorState: "grab"
                  })
              }))
          }
          isEventWithinComponent(e) {
              const relPos = mouse.rel(e);
              return 0 <= relPos[0] && relPos[0] <= this.props.width && 0 <= relPos[1] && relPos[1] <= this.props.height
          }
          componentWillUnmount() {
              this.container.current.removeEventListener("mouseenter", this.onMouseEnter),
              this.container.current.removeEventListener("mouseleave", this.onMouseLeave),
              this.container.current.removeEventListener("click", this.onClick),
              this.container.current.removeEventListener("dblclick", this.onDoubleClick),
              this.container.current.removeEventListener("touchstart", this.onTouchStart),
              this.container.current.removeEventListener("touchend", this.onTouchEnd),
              this.container.current.removeEventListener("touchcancel", this.onTouchCancel),
              this.container.current.removeEventListener("touchmove", this.onTouchMove),
              this.props.sequenceDisableDragging || (this.container.current.removeEventListener("mouseup", this.onMouseUp),
              this.container.current.removeEventListener("mousedown", this.onMouseDown),
              this.container.current.removeEventListener("mousemove", this.onMouseMove)),
              this.stopDragPhase()
          }
          render() {
              const style = DraggingComponent_objectSpread({
                  width: this.props.width,
                  height: this.props.height
              }, this.props.style, {
                  cursor: this.state.mouse.cursorState,
                  position: "relative",
                  flexShrink: 0
              })
                , showModBar = this.props.showModBar && this.state.mouse.isMouseWithin
                , canvasStyle = {
                  position: "absolute",
                  left: 0,
                  top: 0
              }
                , otherProps = Object(omit.a)(this.props, [...this.constructor.propKeys, "tileWidth", "tileHeight", "colorScheme", "nrXTiles", "nrYTiles", "dispatch", "sequences", "fullWidth", "fullHeight", "position", "positionDispatch"]);
              return react_default.a.createElement("div", _extends({
                  style: style,
                  ref: this.container
              }, otherProps), showModBar && react_default.a.createElement(components_ModBar, {
                  style: {
                      position: "absolute",
                      right: 0,
                      opacity: .8
                  }
              }, " Plotly Modbar"), react_default.a.createElement("canvas", {
                  style: canvasStyle,
                  ref: this.canvasBuffers[0],
                  width: this.props.width,
                  height: this.props.height
              }, "Your browser does not seem to support HTML5 canvas."), react_default.a.createElement("canvas", {
                  style: canvasStyle,
                  ref: this.canvasBuffers[1],
                  width: this.props.width,
                  height: this.props.height
              }, "Your browser does not seem to support HTML5 canvas."), react_default.a.createElement(Canvas_FakeScroll, {
                  overflow: this.props.overflow,
                  overflowX: this.props.overflowX,
                  overflowY: this.props.overflowY,
                  positionX: this.props.scrollBarPositionX,
                  positionY: this.props.scrollBarPositionY,
                  width: this.props.width,
                  height: this.props.height,
                  fullWidth: this.props.fullWidth,
                  fullHeight: this.props.fullHeight
              }))
          }
      }
      DraggingComponent_DraggingComponent.displayName = "DraggingComponent",
      DraggingComponent_DraggingComponent.defaultProps = {
          overflow: "auto",
          overflowX: "initial",
          overflowY: "initial",
          scrollBarPositionX: "bottom",
          scrollBarPositionY: "right",
          scrollBarWidth: 5,
          showModBar: !0,
          engine: "canvas",
          sequenceDisableDragging: !1
      },
      DraggingComponent_DraggingComponent.propTypes = {
          overflow: prop_types_default.a.oneOf(["hidden", "auto", "scroll"]),
          overflowX: prop_types_default.a.oneOf(["hidden", "auto", "scroll", "initial"]),
          overflowY: prop_types_default.a.oneOf(["hidden", "auto", "scroll", "initial"]),
          width: prop_types_default.a.number.isRequired,
          height: prop_types_default.a.number.isRequired,
          fullHeight: prop_types_default.a.number.isRequired,
          fullWidth: prop_types_default.a.number.isRequired,
          scrollBarPositionX: prop_types_default.a.oneOf(["top", "bottom"]),
          scrollBarPositionY: prop_types_default.a.oneOf(["left", "right"]),
          showModBar: prop_types_default.a.bool,
          engine: prop_types_default.a.string,
          sequenceDisableDragging: prop_types_default.a.bool
      },
      DraggingComponent_DraggingComponent.__docgenInfo = {
          description: 'Provides dragging support in a canvas for sub-classes.\nSub-classes are expected to implement:\n- drawScene\n- onPositionUpdate(oldPos, newPos)\n\nMoreover, a component\'s viewpoint needs to be passed in via its properties:\n\n  <MyDraggingComponent width="200" height="300" />',
          methods: [{
              name: "drawScene",
              docblock: "Called on every movement to rerender the canvas.",
              modifiers: [],
              params: [],
              returns: null,
              description: "Called on every movement to rerender the canvas."
          }, {
              name: "onPositionUpdate",
              docblock: "Called on every position update.",
              modifiers: [],
              params: [],
              returns: null,
              description: "Called on every position update."
          }, {
              name: "onViewpointChange",
              docblock: "Called every time when the component's dimensions change.",
              modifiers: [],
              params: [],
              returns: null,
              description: "Called every time when the component's dimensions change."
          }, {
              name: "swapContexts",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "draw",
              docblock: "Starts a draw operation by essentially:\n- clearing the current context (the hidden canvas)\n- calling `drawScene` to render the current canvas\n- swapping canvas contexts with `swapContexts`",
              modifiers: [],
              params: [],
              returns: null,
              description: "Starts a draw operation by essentially:\n- clearing the current context (the hidden canvas)\n- calling `drawScene` to render the current canvas\n- swapping canvas contexts with `swapContexts`"
          }, {
              name: "onClick",
              docblock: "To be implemented by its childs.",
              modifiers: [],
              params: [{
                  name: "e"
              }],
              returns: null,
              description: "To be implemented by its childs."
          }, {
              name: "onDoubleClick",
              docblock: "To be implemented by its childs.",
              modifiers: [],
              params: [{
                  name: "e"
              }],
              returns: null,
              description: "To be implemented by its childs."
          }, {
              name: "onMouseDown",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "e",
                  type: null
              }],
              returns: null
          }, {
              name: "onMouseMove",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "e",
                  type: null
              }],
              returns: null
          }, {
              name: "onMouseUp",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "onMouseEnter",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "onMouseLeave",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "onTouchStart",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "e",
                  type: null
              }],
              returns: null
          }, {
              name: "onTouchMove",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "e",
                  type: null
              }],
              returns: null
          }, {
              name: "onTouchEnd",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "onTouchCancel",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "startDragPhase",
              docblock: "Called at the start of a drag action.",
              modifiers: [],
              params: [{
                  name: "e"
              }],
              returns: null,
              description: "Called at the start of a drag action."
          }, {
              name: "stopHoverPhase",
              docblock: "Called whenever the mouse leaves the canvas area.",
              modifiers: [],
              params: [],
              returns: null,
              description: "Called whenever the mouse leaves the canvas area."
          }, {
              name: "stopDragPhase",
              docblock: "Called at the end of a drag action.",
              modifiers: [],
              params: [],
              returns: null,
              description: "Called at the end of a drag action."
          }, {
              name: "isEventWithinComponent",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "e",
                  type: null
              }],
              returns: null
          }],
          displayName: "DraggingComponent",
          props: {
              overflow: {
                  defaultValue: {
                      value: '"auto"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"hidden"',
                          computed: !1
                      }, {
                          value: '"auto"',
                          computed: !1
                      }, {
                          value: '"scroll"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: ""
              },
              overflowX: {
                  defaultValue: {
                      value: '"initial"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"hidden"',
                          computed: !1
                      }, {
                          value: '"auto"',
                          computed: !1
                      }, {
                          value: '"scroll"',
                          computed: !1
                      }, {
                          value: '"initial"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: ""
              },
              overflowY: {
                  defaultValue: {
                      value: '"initial"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"hidden"',
                          computed: !1
                      }, {
                          value: '"auto"',
                          computed: !1
                      }, {
                          value: '"scroll"',
                          computed: !1
                      }, {
                          value: '"initial"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: ""
              },
              scrollBarPositionX: {
                  defaultValue: {
                      value: '"bottom"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"top"',
                          computed: !1
                      }, {
                          value: '"bottom"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: ""
              },
              scrollBarPositionY: {
                  defaultValue: {
                      value: '"right"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"left"',
                          computed: !1
                      }, {
                          value: '"right"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: ""
              },
              scrollBarWidth: {
                  defaultValue: {
                      value: "5",
                      computed: !1
                  },
                  required: !1
              },
              showModBar: {
                  defaultValue: {
                      value: "true",
                      computed: !1
                  },
                  type: {
                      name: "bool"
                  },
                  required: !1,
                  description: ""
              },
              engine: {
                  defaultValue: {
                      value: '"canvas"',
                      computed: !1
                  },
                  type: {
                      name: "string"
                  },
                  required: !1,
                  description: ""
              },
              sequenceDisableDragging: {
                  defaultValue: {
                      value: "false",
                      computed: !1
                  },
                  type: {
                      name: "bool"
                  },
                  required: !1,
                  description: ""
              },
              width: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              height: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              fullHeight: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              },
              fullWidth: {
                  type: {
                      name: "number"
                  },
                  required: !0,
                  description: ""
              }
          }
      };
      var Canvas_DraggingComponent = DraggingComponent_DraggingComponent;
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/Canvas/DraggingComponent.js"] = {
          name: "DraggingComponent",
          docgenInfo: DraggingComponent_DraggingComponent.__docgenInfo,
          path: "src/components/Canvas/DraggingComponent.js"
      });
      var CanvasComponent = __webpack_require__(159);
      class CanvasTilingGrid_CanvasTilingGridComponent extends CanvasComponent.a {
          drawTile({row: row, column: column}) {
              const tileWidth = this.props.tileWidth
                , tileHeight = this.props.tileHeight
                , yPos = tileHeight * (row - this.props.startYTile)
                , xPos = tileWidth * (column - this.props.startXTile);
              if (row >= this.props.sequences.raw.length)
                  return;
              const sequence = this.props.sequences.raw[row].sequence;
              if (column >= sequence.length)
                  return;
              const text = sequence[column];
              if (void 0 !== text) {
                  const colorScheme = this.props.colorScheme.getColor(text, column)
                    , overlayFactor = this.getOverlayFactor(text, column)
                    , key = "".concat(text, "-").concat(colorScheme, "-").concat(overlayFactor)
                    , canvasTile = this.props.residueTileCache.createTile({
                      key: key,
                      tileWidth: tileWidth,
                      tileHeight: tileHeight,
                      create: ({canvas: canvas})=>this.drawResidue({
                          text: text,
                          canvas: canvas,
                          row: row,
                          column: column,
                          colorScheme: colorScheme,
                          overlayFactor: overlayFactor
                      })
                  });
                  this.props.ctx.drawImage(canvasTile, 0, 0, tileWidth, tileHeight, xPos, yPos, tileWidth, tileHeight)
              }
          }
          getOverlayFactor(text, column) {
              if (!this.props.overlayConservation || !this.props.conservation)
                  return 1;
              const raw = this.props.conservation.map[column][text] || 0;
              return Math.floor(5 * raw) / 5
          }
          drawResidue({row: row, column: column, canvas: canvas, colorScheme: colorScheme, text: text, overlayFactor: overlayFactor=1}) {
              if (canvas.globalAlpha = .7 * overlayFactor,
              canvas.fillStyle = colorScheme,
              canvas.fillRect(0, 0, this.props.tileWidth, this.props.tileHeight),
              this.props.tileWidth < 4)
                  return;
              this.props.border && (canvas.globalAlpha = 1,
              canvas.lineWidth = this.props.borderWidth,
              canvas.strokeStyle = this.props.borderStyle,
              canvas.strokeRect(0, 0, this.props.tileWidth, this.props.tileHeight));
              canvas.globalAlpha = Math.min(1, 1 / 6 * this.props.tileWidth + -1 / 6 * 4),
              canvas.fillStyle = this.props.textColor,
              canvas.font = this.props.textFont,
              canvas.textBaseline = "middle",
              canvas.textAlign = "center",
              canvas.fillText(text, this.props.tileWidth / 2, this.props.tileHeight / 2 + 1, this.props.tileWidth)
          }
          draw(props) {
              this.props = props;
              for (let i = this.props.startYTile; i < this.props.endYTile; i++)
                  for (let j = this.props.startXTile; j < this.props.endXTile; j++)
                      this.drawTile({
                          row: i,
                          column: j
                      })
          }
      }
      var CanvasTilingGrid = CanvasTilingGrid_CanvasTilingGridComponent
        , shallowEqual = __webpack_require__(93)
        , minBy = __webpack_require__(928);
      var Canvas_CanvasCache = class CanvasCache_CanvasCache {
          constructor({maxElements: maxElements}={}) {
              this.maxElements = maxElements || 200,
              this.invalidate()
          }
          createTile({key: key, tileWidth: tileWidth, tileHeight: tileHeight, create: create}) {
              if (key in this.cache)
                  return this.cache[key].value;
              if (this.cachedElements >= this.maxElements) {
                  const oldestKey = Object(minBy.a)(Object.keys(this.cache), k=>this.cache[k].insertionTime);
                  delete this.cache[oldestKey]
              }
              const canvas = document.createElement("canvas");
              this.cache[key] = {
                  value: canvas,
                  insertionTime: Date.now()
              },
              canvas.width = tileWidth,
              canvas.height = tileHeight;
              const ctx = canvas.getContext("2d");
              return this.cachedElements++,
              create({
                  canvas: ctx
              }),
              canvas
          }
          updateTileSpecs(spec) {
              return !Object(shallowEqual.a)(spec, this.spec) && (this.invalidate(),
              this.spec = spec,
              !0)
          }
          invalidate() {
              this.cache = {},
              this.spec = {},
              this.cachedElements = 0
          }
      }
        , connect = __webpack_require__(16);
      function roundMod(number, mod) {
          return number - number % mod
      }
      var debug = __webpack_require__(158);
      function _slicedToArray(arr, i) {
          return function _arrayWithHoles(arr) {
              if (Array.isArray(arr))
                  return arr
          }(arr) || function _iterableToArrayLimit(arr, i) {
              if ("undefined" == typeof Symbol || !(Symbol.iterator in Object(arr)))
                  return;
              var _arr = []
                , _n = !0
                , _d = !1
                , _e = void 0;
              try {
                  for (var _s, _i = arr[Symbol.iterator](); !(_n = (_s = _i.next()).done) && (_arr.push(_s.value),
                  !i || _arr.length !== i); _n = !0)
                      ;
              } catch (err) {
                  _d = !0,
                  _e = err
              } finally {
                  try {
                      _n || null == _i.return || _i.return()
                  } finally {
                      if (_d)
                          throw _e
                  }
              }
              return _arr
          }(arr, i) || _unsupportedIterableToArray(arr, i) || function _nonIterableRest() {
              throw new TypeError("Invalid attempt to destructure non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.")
          }()
      }
      function _unsupportedIterableToArray(o, minLen) {
          if (o) {
              if ("string" == typeof o)
                  return _arrayLikeToArray(o, minLen);
              var n = Object.prototype.toString.call(o).slice(8, -1);
              return "Object" === n && o.constructor && (n = o.constructor.name),
              "Map" === n || "Set" === n ? Array.from(n) : "Arguments" === n || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n) ? _arrayLikeToArray(o, minLen) : void 0
          }
      }
      function _arrayLikeToArray(arr, len) {
          (null == len || len > arr.length) && (len = arr.length);
          for (var i = 0, arr2 = Array(len); i < len; i++)
              arr2[i] = arr[i];
          return arr2
      }
      function SequenceViewer_ownKeys(object, enumerableOnly) {
          var keys = Object.keys(object);
          if (Object.getOwnPropertySymbols) {
              var symbols = Object.getOwnPropertySymbols(object);
              enumerableOnly && (symbols = symbols.filter((function(sym) {
                  return Object.getOwnPropertyDescriptor(object, sym).enumerable
              }
              ))),
              keys.push.apply(keys, symbols)
          }
          return keys
      }
      function SequenceViewer_defineProperty(obj, key, value) {
          return key in obj ? Object.defineProperty(obj, key, {
              value: value,
              enumerable: !0,
              configurable: !0,
              writable: !0
          }) : obj[key] = value,
          obj
      }
      class SequenceViewer_SequenceViewerComponent extends Canvas_DraggingComponent {
          constructor(props) {
              super(props),
              SequenceViewer_defineProperty(this, "renderTile", ({row: row, column: column})=>this.tileCache.createTile({
                  key: row + "-" + column,
                  tileWidth: this.props.tileWidth * this.props.xGridSize,
                  tileHeight: this.props.tileHeight * this.props.yGridSize,
                  create: ({canvas: canvas})=>{
                      debug.a && this.redrawnTiles++,
                      this.tilingGridManager.draw(function SequenceViewer_objectSpread(target) {
                          for (var source, i = 1; i < arguments.length; i++)
                              source = null != arguments[i] ? arguments[i] : {},
                              i % 2 ? SequenceViewer_ownKeys(Object(source), !0).forEach((function(key) {
                                  SequenceViewer_defineProperty(target, key, source[key])
                              }
                              )) : Object.getOwnPropertyDescriptors ? Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)) : SequenceViewer_ownKeys(Object(source)).forEach((function(key) {
                                  Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key))
                              }
                              ));
                          return target
                      }({
                          ctx: canvas,
                          startYTile: row,
                          startXTile: column,
                          residueTileCache: this.residueTileCache,
                          endYTile: row + this.props.yGridSize,
                          endXTile: column + this.props.xGridSize
                      }, Object(pick.a)(this.props, ["sequences", "colorScheme", "textFont", "textColor", "tileHeight", "tileWidth", "border", "borderWidth", "borderColor", "overlayConservation", "conservation"])))
                  }
              })),
              SequenceViewer_defineProperty(this, "onPositionUpdate", (oldPos,newPos)=>{
                  const relativeMovement = {
                      xMovement: oldPos[0] - newPos[0],
                      yMovement: oldPos[1] - newPos[1]
                  };
                  this.sendEvent("onPositionUpdate", newPos),
                  this.props.positionDispatch(Object(positionReducers.b)(relativeMovement))
              }
              ),
              SequenceViewer_defineProperty(this, "updateScrollPosition", ()=>{
                  this.draw()
              }
              ),
              SequenceViewer_defineProperty(this, "onMouseMove", e=>{
                  if (void 0 === this.isInDragPhase && this.hasOnMouseMoveProps()) {
                      const eventData = this.currentPointerPosition(e)
                        , lastValue = this.currentMouseSequencePosition;
                      if (!Object(isEqual.a)(lastValue, eventData) && (void 0 !== lastValue && this.sendEvent("onResidueMouseLeave", lastValue),
                      this.currentMouseSequencePosition = eventData,
                      this.sendEvent("onResidueMouseEnter", eventData),
                      this.props.features && 0 < this.props.features.length)) {
                          const lastMouseOverFeatureIds = this.mouseOverFeatureIds || []
                            , mouseOverFeatureIds = this.sequencePositionToFeatureIds(eventData);
                          Object(isEqual.a)(lastMouseOverFeatureIds, mouseOverFeatureIds) || (this.mouseOverFeatureIds = mouseOverFeatureIds,
                          super.draw())
                      }
                  }
                  super.onMouseMove(e)
              }
              ),
              SequenceViewer_defineProperty(this, "onMouseLeave", e=>{
                  this.sendEvent("onResidueMouseLeave", this.currentMouseSequencePosition),
                  this.currentMouseSequencePosition = void 0,
                  this.mouseOverFeatureIds && (this.mouseOverFeatureIds = void 0,
                  super.draw()),
                  super.onMouseLeave(e)
              }
              ),
              SequenceViewer_defineProperty(this, "onClick", e=>{
                  if (!this.mouseHasMoved) {
                      const eventData = this.currentPointerPosition(e);
                      this.sendEvent("onResidueClick", eventData),
                      this.props.features && this.sequencePositionToFeatureIds(eventData).forEach(id=>{
                          this.sendEvent("onFeatureClick", {
                              event: e,
                              id: id
                          })
                      }
                      )
                  }
                  super.onClick(e)
              }
              ),
              SequenceViewer_defineProperty(this, "onDoubleClick", e=>{
                  const eventData = this.currentPointerPosition(e);
                  this.sendEvent("onResidueDoubleClick", eventData),
                  super.onDoubleClick(e)
              }
              ),
              this.tileCache = new Canvas_CanvasCache,
              this.residueTileCache = new Canvas_CanvasCache,
              this.tilingGridManager = new CanvasTilingGrid
          }
          hasOnMouseMoveProps() {
              return void 0 !== this.props.onResidueMouseEnter || void 0 !== this.props.onResidueMouseLeave || this.props.features && 0 < this.props.features.length
          }
          componentDidMount() {
              this.hasOnMouseMoveProps() && this.container.current.addEventListener("mousemove", this.onMouseMove),
              super.componentDidMount()
          }
          componentWillUnmount() {
              this.hasOnMouseMoveProps() && this.container.current.removeEventListener("mousemove", this.onMouseMove),
              super.componentWillUnmount()
          }
          drawScene() {
              const positions = this.getTilePositions();
              if (this.updateTileSpecs(),
              debug.a && (this.redrawStarted = Date.now(),
              this.redrawnTiles = 0),
              this.drawTiles(positions),
              this.drawHighlightedRegions(),
              this.ctx && this.ctx.canvas.dispatchEvent(new CustomEvent("drawCompleted",{
                  bubbles: !0
              })),
              debug.a) {
                  const elapsed = Date.now() - this.redrawStarted;
                  5 < elapsed && console.log("Took ".concat(elapsed, " msecs to redraw for ").concat(positions.startXTile, " ").concat(positions.startYTile, " (redrawnTiles: ").concat(this.redrawnTiles, ")"))
              }
          }
          getTilePositions() {
              const startXTile = Math.max(0, this.props.position.currentViewSequencePosition - this.props.cacheElements)
                , startYTile = Math.max(0, this.props.position.currentViewSequence - this.props.cacheElements)
                , endYTile = Math.min(this.props.sequences.length, startYTile + this.props.nrYTiles + 2 * this.props.cacheElements);
              return {
                  startXTile: startXTile,
                  startYTile: startYTile,
                  endXTile: Math.min(this.props.sequences.maxLength, startXTile + this.props.nrXTiles + 2 * this.props.cacheElements),
                  endYTile: endYTile
              }
          }
          drawTiles({startXTile: startXTile, startYTile: startYTile, endXTile: endXTile, endYTile: endYTile}) {
              const xGridSize = this.props.xGridSize
                , yGridSize = this.props.yGridSize
                , startY = roundMod(startYTile, yGridSize)
                , startX = roundMod(startXTile, xGridSize);
              for (let i = startY; i < endYTile; i += yGridSize)
                  for (let j = startX; j < endXTile; j += xGridSize) {
                      const canvas = this.renderTile({
                          row: i,
                          column: j,
                          canvas: this.ctx
                      })
                        , width = xGridSize * this.props.tileWidth
                        , height = yGridSize * this.props.tileHeight
                        , yPos = (i - this.props.position.currentViewSequence) * this.props.tileHeight + this.props.position.yPosOffset
                        , xPos = (j - this.props.position.currentViewSequencePosition) * this.props.tileWidth + this.props.position.xPosOffset;
                      this.ctx.drawImage(canvas, 0, 0, width, height, xPos, yPos, width, height)
                  }
          }
          drawHighlightedRegions() {
              if (this.props.highlight)
                  if (Array.isArray(this.props.highlight)) {
                      var _step, _iterator = function _createForOfIteratorHelper(o) {
                          if ("undefined" == typeof Symbol || null == o[Symbol.iterator]) {
                              if (Array.isArray(o) || (o = _unsupportedIterableToArray(o))) {
                                  var i = 0
                                    , F = function() {};
                                  return {
                                      s: F,
                                      n: function n() {
                                          return i >= o.length ? {
                                              done: !0
                                          } : {
                                              done: !1,
                                              value: o[i++]
                                          }
                                      },
                                      e: function e(_e2) {
                                          throw _e2
                                      },
                                      f: F
                                  }
                              }
                              throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.")
                          }
                          var it, err, normalCompletion = !0, didErr = !1;
                          return {
                              s: function s() {
                                  it = o[Symbol.iterator]()
                              },
                              n: function n() {
                                  var step = it.next();
                                  return normalCompletion = step.done,
                                  step
                              },
                              e: function e(_e3) {
                                  didErr = !0,
                                  err = _e3
                              },
                              f: function f() {
                                  try {
                                      normalCompletion || null == it.return || it.return()
                                  } finally {
                                      if (didErr)
                                          throw err
                                  }
                              }
                          }
                      }(this.props.highlight);
                      try {
                          for (_iterator.s(); !(_step = _iterator.n()).done; ) {
                              const h = _step.value;
                              this.drawHighligtedRegion(h)
                          }
                      } catch (err) {
                          _iterator.e(err)
                      } finally {
                          _iterator.f()
                      }
                  } else
                      this.drawHighligtedRegion(this.props.highlight);
              this.props.features && this.props.features.forEach(feature=>{
                  this.drawHighligtedRegion(feature)
              }
              )
          }
          drawHighligtedRegion(region) {
              var _this$mouseOverFeatur;
              if (!this.ctx || !region)
                  return;
              const regionWidth = this.props.tileWidth * (1 + region.residues.to - region.residues.from)
                , regionHeight = this.props.tileHeight * (1 + region.sequences.to - region.sequences.from)
                , yPosFrom = (region.sequences.from - this.props.position.currentViewSequence) * this.props.tileHeight + this.props.position.yPosOffset
                , xPosFrom = (region.residues.from - 1 - this.props.position.currentViewSequencePosition) * this.props.tileWidth + this.props.position.xPosOffset
                , canvas = document.createElement("canvas");
              canvas.width = regionWidth,
              canvas.height = regionHeight;
              const ctx = canvas.getContext("2d")
                , mouseOver = null === (_this$mouseOverFeatur = this.mouseOverFeatureIds) || void 0 === _this$mouseOverFeatur ? void 0 : _this$mouseOverFeatur.some(id=>id === region.id);
              ctx.globalAlpha = .3,
              ctx.fillStyle = mouseOver ? region.mouseOverFillColor || "green" : region.fillColor || "#9999FF",
              ctx.fillRect(0, 0, regionWidth, regionHeight),
              ctx.globalAlpha = 1,
              ctx.strokeStyle = mouseOver ? region.mouseOverBorderColor || "black " : region.borderColor || "777700",
              ctx.lineWidth = "4",
              ctx.rect(0, 0, regionWidth, regionHeight),
              ctx.stroke(),
              this.ctx.drawImage(canvas, 0, 0, regionWidth, regionHeight, xPosFrom, yPosFrom, regionWidth, regionHeight)
          }
          positionToSequence(pos) {
              const sequences = this.props.sequences.raw
                , seqNr = Object(clamp.a)(Object(floor.a)((this.props.position.yPos + pos.yPos) / this.props.tileHeight), 0, sequences.length - 1)
                , sequence = sequences[seqNr]
                , position = Object(clamp.a)(Object(floor.a)((this.props.position.xPos + pos.xPos) / this.props.tileWidth), 0, sequence.sequence.length - 1);
              return {
                  i: seqNr,
                  sequence: sequence,
                  position: position,
                  residue: sequence.sequence[position]
              }
          }
          sequencePositionToFeatureIds(sequencePosition) {
              if (this.props.features)
                  return this.props.features.filter(feature=>feature.id && sequencePosition.position >= feature.residues.from - 1 && sequencePosition.position <= feature.residues.to - 1 && sequencePosition.i >= feature.sequences.from && sequencePosition.i <= feature.sequences.to).map(feature=>feature.id)
          }
          currentPointerPosition(e) {
              const _Mouse$rel2 = _slicedToArray(mouse.rel(e), 2)
                , x = _Mouse$rel2[0]
                , y = _Mouse$rel2[1];
              return this.positionToSequence({
                  xPos: x,
                  yPos: y
              })
          }
          sendEvent(name, data) {
              void 0 !== this.props[name] && this.props[name](data)
          }
          componentWillUnmount() {
              this.tileCache.invalidate(),
              this.residueTileCache.invalidate()
          }
          updateTileSpecs() {
              const tileAttributes = ["tileWidth", "tileHeight", "colorScheme", "textFont", "borderColor", "overlayConservation", "conservation"];
              this.tileCache.updateTileSpecs(Object(pick.a)(this.props, [...tileAttributes, "xGridSize", "yGridSize", "sequences"])),
              this.residueTileCache.updateTileSpecs(Object(pick.a)(this.props, tileAttributes))
          }
          render() {
              return super.render()
          }
      }
      SequenceViewer_SequenceViewerComponent.displayName = "SequenceViewerComponent",
      SequenceViewer_SequenceViewerComponent.defaultProps = {
          showModBar: !1,
          xGridSize: 10,
          yGridSize: 10,
          border: !1,
          borderColor: "black",
          borderWidth: 1,
          cacheElements: 20,
          textColor: "black",
          textFont: "18px Arial",
          overflow: "hidden",
          overflowX: "auto",
          overflowY: "auto",
          scrollBarPositionX: "bottom",
          scrollBarPositionY: "right",
          overlayConservation: !1,
          conservation: null,
          sequenceDisableDragging: !1,
          highlight: null
      },
      SequenceViewer_SequenceViewerComponent.propTypes = {
          showModBar: prop_types_default.a.bool,
          onResidueMouseEnter: prop_types_default.a.func,
          onResidueMouseLeave: prop_types_default.a.func,
          onResidueClick: prop_types_default.a.func,
          onFeatureClick: prop_types_default.a.func,
          highlight: prop_types_default.a.oneOfType([prop_types_default.a.object, prop_types_default.a.array]),
          features: prop_types_default.a.arrayOf(prop_types_default.a.object),
          onResidueDoubleClick: prop_types_default.a.func,
          xGridSize: prop_types_default.a.number.isRequired,
          yGridSize: prop_types_default.a.number.isRequired,
          cacheElements: prop_types_default.a.number.isRequired,
          border: prop_types_default.a.bool,
          borderColor: prop_types_default.a.string,
          borderWidth: prop_types_default.a.number,
          textColor: prop_types_default.a.string,
          textFont: prop_types_default.a.string,
          overflow: prop_types_default.a.oneOf(["hidden", "auto", "scroll"]),
          overflowX: prop_types_default.a.oneOf(["hidden", "auto", "scroll", "initial"]),
          overflowY: prop_types_default.a.oneOf(["hidden", "auto", "scroll", "initial"]),
          scrollBarPositionX: prop_types_default.a.oneOf(["top", "bottom"]),
          scrollBarPositionY: prop_types_default.a.oneOf(["left", "right"]),
          overlayConservation: prop_types_default.a.bool,
          onPositionUpdate: prop_types_default.a.func,
          sequenceDisableDragging: prop_types_default.a.bool
      },
      SequenceViewer_SequenceViewerComponent.propKeys = Object.keys(SequenceViewer_SequenceViewerComponent.propTypes);
      const SV = Object(withPositionStore.a)(SequenceViewer_SequenceViewerComponent, {
          withX: !0,
          withY: !0
      });
      __webpack_exports__.a = Object(connect.a)(state=>{
          const width = Math.min(state.props.width, state.sequences.maxLength * state.props.tileWidth)
            , height = Math.min(state.props.height, state.sequences.length * state.props.tileHeight);
          return {
              sequences: state.sequences,
              width: width,
              height: height,
              highlight: state.props.highlight,
              tileWidth: state.props.tileWidth,
              tileHeight: state.props.tileHeight,
              colorScheme: state.props.colorScheme,
              overlayConservation: state.props.overlayConservation,
              conservation: state.props.overlayConservation ? state.props.conservation : null,
              nrXTiles: state.sequenceStats.nrXTiles,
              nrYTiles: state.sequenceStats.nrYTiles,
              fullWidth: state.sequenceStats.fullWidth,
              fullHeight: state.sequenceStats.fullHeight
          }
      }
      )(SV);
      SequenceViewer_SequenceViewerComponent.__docgenInfo = {
          description: "Component to draw the main sequence alignment.",
          methods: [{
              name: "hasOnMouseMoveProps",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "drawScene",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "getTilePositions",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "renderTile",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "{ row, column }",
                  type: null
              }],
              returns: null
          }, {
              name: "drawTiles",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "{ startXTile, startYTile, endXTile, endYTile }",
                  type: null
              }],
              returns: null
          }, {
              name: "drawHighlightedRegions",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "drawHighligtedRegion",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "region",
                  type: null
              }],
              returns: null
          }, {
              name: "onPositionUpdate",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "oldPos",
                  type: null
              }, {
                  name: "newPos",
                  type: null
              }],
              returns: null
          }, {
              name: "positionToSequence",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "pos",
                  type: null
              }],
              returns: null
          }, {
              name: "sequencePositionToFeatureIds",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "sequencePosition",
                  type: null
              }],
              returns: null
          }, {
              name: "updateScrollPosition",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }, {
              name: "currentPointerPosition",
              docblock: "Returns the position of the mouse position relative to the sequences",
              modifiers: [],
              params: [{
                  name: "e"
              }],
              returns: null,
              description: "Returns the position of the mouse position relative to the sequences"
          }, {
              name: "sendEvent",
              docblock: "Only sends an event if the actual function is set.",
              modifiers: [],
              params: [{
                  name: "name"
              }, {
                  name: "data"
              }],
              returns: null,
              description: "Only sends an event if the actual function is set."
          }, {
              name: "onMouseMove",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "e",
                  type: null
              }],
              returns: null
          }, {
              name: "onMouseLeave",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "e",
                  type: null
              }],
              returns: null
          }, {
              name: "onClick",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "e",
                  type: null
              }],
              returns: null
          }, {
              name: "onDoubleClick",
              docblock: null,
              modifiers: [],
              params: [{
                  name: "e",
                  type: null
              }],
              returns: null
          }, {
              name: "updateTileSpecs",
              docblock: null,
              modifiers: [],
              params: [],
              returns: null
          }],
          displayName: "SequenceViewerComponent",
          props: {
              showModBar: {
                  defaultValue: {
                      value: "false",
                      computed: !1
                  },
                  type: {
                      name: "bool"
                  },
                  required: !1,
                  description: "Show the custom ModBar"
              },
              xGridSize: {
                  defaultValue: {
                      value: "10",
                      computed: !1
                  },
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Number of residues to cluster in one tile (x-axis) (default: 10)"
              },
              yGridSize: {
                  defaultValue: {
                      value: "10",
                      computed: !1
                  },
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Number of residues to cluster in one tile (y-axis) (default: 10)"
              },
              border: {
                  defaultValue: {
                      value: "false",
                      computed: !1
                  },
                  type: {
                      name: "bool"
                  },
                  required: !1,
                  description: "Whether to draw a border."
              },
              borderColor: {
                  defaultValue: {
                      value: '"black"',
                      computed: !1
                  },
                  type: {
                      name: "string"
                  },
                  required: !1,
                  description: "Color of the border. Name, hex or RGB value."
              },
              borderWidth: {
                  defaultValue: {
                      value: "1",
                      computed: !1
                  },
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Width of the border."
              },
              cacheElements: {
                  defaultValue: {
                      value: "20",
                      computed: !1
                  },
                  type: {
                      name: "number"
                  },
                  required: !1,
                  description: "Number of residues to prerender outside of the visible viewbox."
              },
              textColor: {
                  defaultValue: {
                      value: '"black"',
                      computed: !1
                  },
                  type: {
                      name: "string"
                  },
                  required: !1,
                  description: "Color of the text residue letters (name, hex or RGB value)"
              },
              textFont: {
                  defaultValue: {
                      value: '"18px Arial"',
                      computed: !1
                  },
                  type: {
                      name: "string"
                  },
                  required: !1,
                  description: "Font to use when drawing the individual residues."
              },
              overflow: {
                  defaultValue: {
                      value: '"hidden"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"hidden"',
                          computed: !1
                      }, {
                          value: '"auto"',
                          computed: !1
                      }, {
                          value: '"scroll"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: "What should happen if content overflows."
              },
              overflowX: {
                  defaultValue: {
                      value: '"auto"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"hidden"',
                          computed: !1
                      }, {
                          value: '"auto"',
                          computed: !1
                      }, {
                          value: '"scroll"',
                          computed: !1
                      }, {
                          value: '"initial"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: 'What should happen if x-axis content overflows (overwrites "overflow")'
              },
              overflowY: {
                  defaultValue: {
                      value: '"auto"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"hidden"',
                          computed: !1
                      }, {
                          value: '"auto"',
                          computed: !1
                      }, {
                          value: '"scroll"',
                          computed: !1
                      }, {
                          value: '"initial"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: 'What should happen if y-axis content overflows (overwrites "overflow")'
              },
              scrollBarPositionX: {
                  defaultValue: {
                      value: '"bottom"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"top"',
                          computed: !1
                      }, {
                          value: '"bottom"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: 'X Position of the scroll bar ("top or "bottom")'
              },
              scrollBarPositionY: {
                  defaultValue: {
                      value: '"right"',
                      computed: !1
                  },
                  type: {
                      name: "enum",
                      value: [{
                          value: '"left"',
                          computed: !1
                      }, {
                          value: '"right"',
                          computed: !1
                      }]
                  },
                  required: !1,
                  description: 'Y Position of the scroll bar ("left" or "right")'
              },
              overlayConservation: {
                  defaultValue: {
                      value: "false",
                      computed: !1
                  },
                  type: {
                      name: "bool"
                  },
                  required: !1,
                  description: "The conservation data can be used to define an overlay that\ndefines the opacity of the background color of each residue."
              },
              conservation: {
                  defaultValue: {
                      value: "null",
                      computed: !1
                  },
                  required: !1
              },
              sequenceDisableDragging: {
                  defaultValue: {
                      value: "false",
                      computed: !1
                  },
                  type: {
                      name: "bool"
                  },
                  required: !1,
                  description: ""
              },
              highlight: {
                  defaultValue: {
                      value: "null",
                      computed: !1
                  },
                  type: {
                      name: "union",
                      value: [{
                          name: "object"
                      }, {
                          name: "array"
                      }]
                  },
                  required: !1,
                  description: "Displays a highlight"
              },
              onResidueMouseEnter: {
                  type: {
                      name: "func"
                  },
                  required: !1,
                  description: "Callback fired when the mouse pointer is entering a residue."
              },
              onResidueMouseLeave: {
                  type: {
                      name: "func"
                  },
                  required: !1,
                  description: "Callback fired when the mouse pointer is leaving a residue."
              },
              onResidueClick: {
                  type: {
                      name: "func"
                  },
                  required: !1,
                  description: "Callback fired when the mouse pointer clicked a residue."
              },
              onFeatureClick: {
                  type: {
                      name: "func"
                  },
                  required: !1,
                  description: "Callback fired when the mouse pointer clicked a feature."
              },
              features: {
                  type: {
                      name: "arrayOf",
                      value: {
                          name: "object"
                      }
                  },
                  required: !1,
                  description: "An array of features which can be clicked"
              },
              onResidueDoubleClick: {
                  type: {
                      name: "func"
                  },
                  required: !1,
                  description: "Callback fired when the mouse pointer clicked a residue."
              },
              onPositionUpdate: {
                  type: {
                      name: "func"
                  },
                  required: !1,
                  description: ""
              }
          }
      },
      "undefined" != typeof STORYBOOK_REACT_CLASSES && (STORYBOOK_REACT_CLASSES["src/components/Canvas/SequenceViewer.js"] = {
          name: "SequenceViewerComponent",
          docgenInfo: SequenceViewer_SequenceViewerComponent.__docgenInfo,
          path: "src/components/Canvas/SequenceViewer.js"
      })
  }
}, [[470, 1, 2]]]);
//# sourceMappingURL=main.a72dfd6910bd0b973489.bundle.js.map
