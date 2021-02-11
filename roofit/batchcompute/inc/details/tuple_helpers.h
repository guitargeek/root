#include<tuple>

namespace RooBatchCompute {

  namespace tuple_helpers {

    template <class... Args, std::size_t... Is>
    constexpr auto tuple_init_helper(std::tuple<Args...> tp, std::index_sequence<Is...>)
    {
        return std::make_tuple(std::get<Is>(tp)...);
    }

    template <class... Args>
    constexpr auto tuple_init(std::tuple<Args...> tp)
    {
        return tuple_init_helper(tp, std::make_index_sequence<sizeof...(Args) - 1>{});
    }

    template <class... Args, std::size_t... Is>
    constexpr auto tuple_tail_helper(std::tuple<Args...> tp, std::index_sequence<Is...>)
    {
        return std::make_tuple(std::get<Is+1>(tp)...);
    }

    template <class... Args>
    constexpr auto tuple_tail(std::tuple<Args...> tp)
    {
        return tuple_tail_helper(tp, std::make_index_sequence<sizeof...(Args) - 1>{});
    }

    static_assert(tuple_init(std::make_tuple(1))       == std::make_tuple());
    static_assert(tuple_init(std::make_tuple(1, 2))    == std::make_tuple(1));
    static_assert(tuple_init(std::make_tuple(1, 2, 3)) == std::make_tuple(1, 2));

    static_assert(tuple_tail(std::make_tuple(1))       == std::make_tuple());
    static_assert(tuple_tail(std::make_tuple(1, 2))    == std::make_tuple(2));
    static_assert(tuple_tail(std::make_tuple(1, 2, 3)) == std::make_tuple(2, 3));

  }

}
