#include <atomic>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <optional>
#include <thread>

#include <SFML/Graphics.hpp>

namespace bh
{

float THETA;
float GRAVITY_CONSTANT;
float TIME_STEP;
float SOFTENING;

struct point_t
{
  float mass;
  sf::Vector2f position;
  sf::Vector2f velocity;
};

static inline bh::point_t
point_init (float mass, const sf::Vector2f &position,
            const sf::Vector2f &velocity = { 0.f, 0.f })
{
  return (
      bh::point_t){ .mass = mass, .position = position, .velocity = velocity };
}

struct quad_node_t
{
  alignas (8) float total_mass{ 0.f };
  sf::Vector2f center_of_mass{ 0.f, 0.f };
  sf::FloatRect boundary{};
  std::optional<bh::point_t> point{};
  bh::quad_node_t *children[4]{ 0, 0, 0, 0 };
};

static inline bh::quad_node_t *
quad_node_init (const sf::FloatRect &boundary)
{
  auto *node = new bh::quad_node_t{};
  return node->boundary = boundary, node;
}

static inline void
quad_node_free (bh::quad_node_t *node)
{
  if (node == NULL)
    return;

  for (auto child : node->children)
    quad_node_free (child);

  delete node;
}

static inline bool
quad_node_is_leaf (const bh::quad_node_t &node)
{
  return std::all_of (node.children, node.children + 4,
                      [] (const auto &child) { return child == NULL; });
}

static inline void
quad_node_subdivide (bh::quad_node_t *node)
{
  const sf::Vector2f center{ node->boundary.left + node->boundary.width / 2,
                             node->boundary.top + node->boundary.height / 2 };

  const sf::FloatRect quadrants[4]
      = { { node->boundary.left, node->boundary.top, node->boundary.width / 2,
            node->boundary.height / 2 },
          { center.x, node->boundary.top, node->boundary.width / 2,
            node->boundary.height / 2 },
          { node->boundary.left, center.y, node->boundary.width / 2,
            node->boundary.height / 2 },
          { center.x, center.y, node->boundary.width / 2,
            node->boundary.height / 2 } };

  unsigned char index = 0;
  for (auto &child : node->children)
    child = bh::quad_node_init (quadrants[(index++) % 4]);
}

static inline void
quad_node_insert (bh::quad_node_t *node, const bh::point_t &point)
{
  if (!node->boundary.contains (point.position))
    return;

  if (bh::quad_node_is_leaf (*node))
    {
      if (!node->point.has_value ())
        return node->point = point, (void)0;

      bh::quad_node_subdivide (node);

      const bh::point_t save = node->point.value ();
      node->point.reset ();

      for (auto child : node->children)
        bh::quad_node_insert (child, save);
    }

  for (auto child : node->children)
    bh::quad_node_insert (child, point);
}

static inline void
quad_node_compute_mass (bh::quad_node_t *node)
{
  if (bh::quad_node_is_leaf (*node))
    {
      if (node->point.has_value ())
        {
          node->center_of_mass = node->point->position;
          node->total_mass = node->point->mass;
        }

      return;
    }

  node->center_of_mass = { 0.f, 0.f };
  node->total_mass = 0;

  for (auto child : node->children)
    {
      bh::quad_node_compute_mass (child);
      node->total_mass += child->total_mass;
      node->center_of_mass += child->center_of_mass * child->total_mass;
    }

  if (node->total_mass > 0)
    node->center_of_mass /= node->total_mass;
}

static inline void
quad_node_compute_force (const quad_node_t &node, point_t *point)
{
  if (node.total_mass == 0 || point->position == node.center_of_mass)
    return;

  const sf::Vector2f direction = node.center_of_mass - point->position;
  const float distance
      = std::sqrt (direction.x * direction.x + direction.y * direction.y
                   + bh::SOFTENING * bh::SOFTENING);

  const float ratio = node.boundary.width / distance;
  if (bh::quad_node_is_leaf (node) || ratio < THETA)
    {
      const float force
          = bh::GRAVITY_CONSTANT * node.total_mass * point->mass
            / (distance * distance + bh::SOFTENING * bh::SOFTENING);
      const sf::Vector2f acceleration
          = direction / distance * force / point->mass;
      point->velocity += acceleration * bh::TIME_STEP;
    }
  else
    {
      for (auto child : node.children)
        bh::quad_node_compute_force (*child, point);
    }
}

}

#define QT_SIZE 160000

void
push_galaxy (std::vector<bh::point_t> &points, int n, float inital_radius,
             float speed, float center_x, float center_y,
             float base_velocity_x, float base_velocity_y,
             float mass)
{
  for (int i = 0; i < n; ++i)
    {
      float angle = static_cast<float> (std::rand () % 360) * (M_PI / 180.0f);
      float radius = static_cast<float> (std::rand ()) / RAND_MAX;
      radius = sqrtf (radius) * inital_radius;

      float x = center_x + cosf (angle) * radius;
      float y = center_y + sinf (angle) * radius;

      float dx = center_x - x, dy = center_y - y;
      float normal_angle = atan2f (dy, dx) - M_PI / 2;

      points.emplace_back (bh::point_init (
          mass, { x, y },
          {
              base_velocity_x
                  + cosf (normal_angle) * speed * (radius / inital_radius),
              base_velocity_y
                  + sinf (normal_angle) * speed * (radius / inital_radius),
          }));
    }
}

int
main ()
{
  srand(time(nullptr));
  sf::RenderWindow window{ sf::VideoMode{ 800, 800 }, "Barnes-Hut Simulation",
                           sf::Style::Titlebar,
                           sf::ContextSettings{ 24, 8, 8 } };
  window.setFramerateLimit (60);
  window.setPosition ({ 1920 / 2 - 400, 1080 / 2 - 400 });

  sf::VertexArray vao{ sf::Points };

  std::vector<bh::point_t> points{};

  bh::THETA = 0.5f;
  bh::GRAVITY_CONSTANT = 1.0f;
  bh::TIME_STEP = 1.0f;
  bh::SOFTENING = 1.0f;

  push_galaxy (points, 100'000, 400, 12, 0, 0, 0, 0, 1.0);

  std::mutex points_mutex;
  std::atomic<bool> running{ true };

  std::vector<bh::point_t> points_previous = points;
  std::vector<bh::point_t> points_current = points;

  std::vector<bh::point_t> render_previous = points;
  std::vector<bh::point_t> render_current = points;

  std::atomic<bool> update_done = 0;
  std::atomic<bool> do_update = 1;

  auto last_sim_update = std::chrono::steady_clock::now ();
  auto prev_sim_time = std::chrono::steady_clock::now ();
  float sim_update_interval;

  sf::View view (sf::FloatRect (-400, -400, 800, 800));
  float zoom_level = 1.5;

  bool do_interpolate = false;

  view.zoom (zoom_level);

  std::thread sim_thread ([&] () {
    while (running.load ())
      {
        while (!do_update.load ())
          {
            if (!running.load ())
              return;
            std::this_thread::sleep_for (std::chrono::microseconds (1));
          }

        auto start = std::chrono::steady_clock::now ();
        std::vector<bh::point_t> local_points;
        {
          std::lock_guard<std::mutex> lock (points_mutex);
          local_points = points_current;
        }

        bh::quad_node_t *root = bh::quad_node_init (
            { -QT_SIZE, -QT_SIZE, QT_SIZE * 2, QT_SIZE * 2 });
        for (const auto &p : local_points)
          bh::quad_node_insert (root, p);
        bh::quad_node_compute_mass (root);

#pragma omp parallel for
        for (size_t i = 0; i < local_points.size (); ++i)
          {
            bh::quad_node_compute_force (*root, &local_points[i]);
            local_points[i].position
                += local_points[i].velocity * bh::TIME_STEP;
          }

        bh::quad_node_free (root);

        auto now = std::chrono::steady_clock::now ();

        float delta
            = std::chrono::duration<float> (now - prev_sim_time).count ();
        prev_sim_time = now;

        {
          std::lock_guard<std::mutex> lock (points_mutex);

          std::swap (points_previous, points_current);
          std::swap (points_current, local_points);

          last_sim_update = now;
          sim_update_interval = delta;
        }

        auto end = std::chrono::steady_clock::now ();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
            end - start);

        update_done.store (1);
        printf ("\tupdate %ldms\n", duration.count ());
      }
  });

  sf::Clock delta, fps_clock;
  float fps_avg = 0;
  int frame = 0;

  while (window.isOpen ())
    {
    const float dt = delta.restart ().asSeconds ();

    sf::Event event{};
    while (window.pollEvent (event))
      {
        if (event.type == sf::Event::Closed)
          window.close ();

        if (event.type == sf::Event::MouseWheelScrolled)
          {
            if (event.mouseWheelScroll.delta > 0)
              {
                zoom_level /= 1.1f;
                view.zoom (1.f / 1.1f);
              }
            else if (event.mouseWheelScroll.delta < 0)
              {
                zoom_level *= 1.1f;
                view.zoom (1.1f);
              }
          }

        if (event.type == sf::Event::KeyPressed)
          {
            if (event.key.code == sf::Keyboard::Tab)
              {
                do_interpolate = !do_interpolate;
                printf ("do_interpolate=%d\n", do_interpolate);
              }
          }
      }

    const float pan_speed = zoom_level * 300.f * dt;

    if (sf::Keyboard::isKeyPressed (sf::Keyboard::W))
      view.move (0, -pan_speed);
    if (sf::Keyboard::isKeyPressed (sf::Keyboard::S))
      view.move (0, pan_speed);
    if (sf::Keyboard::isKeyPressed (sf::Keyboard::A))
      view.move (-pan_speed, 0);
    if (sf::Keyboard::isKeyPressed (sf::Keyboard::D))
      view.move (pan_speed, 0);

    window.setView (view);
    window.clear ({ 10, 10, 10 });

    auto now = std::chrono::steady_clock::now ();

    if (update_done.load ())
      {
        auto start = std::chrono::steady_clock::now ();

        do_update.store (0);

        render_previous = points_previous;
        render_current = points_current;
        update_done.store (0);

        auto end = std::chrono::steady_clock::now ();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
            end - start);

        printf ("\tcopy   %ldms\n", duration.count ());

        do_update.store (1);
      }

    float alpha = 0.f;
    {
      auto elapsed
          = std::chrono::duration<float> (now - last_sim_update).count ();
      if (do_interpolate)
        alpha = std::clamp (elapsed / sim_update_interval, 0.f, 1.f);

      vao.clear ();

      for (size_t i = 0; i < render_current.size (); ++i)
        {
          const auto &prev = render_previous[i].position;
          const auto &curr = render_current[i].position;

          sf::Vector2f interp_pos = prev + (curr - prev) * alpha;

          sf::Uint8 r = static_cast<sf::Uint8> (92);
          sf::Uint8 g = static_cast<sf::Uint8> (106);
          sf::Uint8 b = static_cast<sf::Uint8> (114);

          sf::Uint8 a = static_cast<sf::Uint8> (128);

          vao.append (sf::Vertex (interp_pos, sf::Color (r, g, b, a)));
        }
    }

    window.draw (vao, sf::RenderStates (sf::BlendAdd));

    sf::RectangleShape shape;
    shape.setSize ({ QT_SIZE * 2, QT_SIZE * 2 });
    shape.setPosition ({ -QT_SIZE, -QT_SIZE });
    shape.setFillColor (sf::Color::Transparent);
    shape.setOutlineColor (sf::Color::White);
    shape.setOutlineThickness (zoom_level);
    window.draw (shape);

    window.display ();

    fps_avg += (frame++ == 0) ? 0 : 1.0f / dt;

    if (fps_clock.getElapsedTime ().asSeconds () > 1.0)
      {
        printf ("%.2f (%.2f)\n", 1.0 / dt, fps_avg / frame);
        fps_clock.restart ();
      }
    }

  running = false;
  sim_thread.join ();

  return 0;
}

